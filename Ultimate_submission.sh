#!/bin/bash

# ######################### #
#     Start of code         #
# ######################### #

# Dictionaries
declare -a errorcodes=( input/input_parsing.F:187 input/input_parsing.F:248 )


# Flags 
corrected='false'
while getopts 'c' flag
do
    case "${flag}" in
        c) corrected='true' ;; # if restarting calculation after fixing a job ending error, 
                               # use the flag -c in your sbatch submission command.
    esac
done

# Variables 
File="LDOS_PBE0_GEO_OPT.log" # log file name 
inp="interface.inp" # inp input file name 
restart="SiC-SiO2_scheme_GEO_OPT-1.restart" # restart input file name
dir="inter_calc_start_from_files" # intermediate input/restart files directory name
xyz_="SiC-SiO2_scheme_GEO_OPT-pos-1.xyz"

atomic_SCF="SCF_GUESS  ATOMIC"
restart_SCF="SCF_GUESS  RESTART"
WFNname="WFN_RESTART_FILE_NAME"
COORDformat="COORD_FILE_FORMAT"
COORDfile="COORD_FILE_NAME"

# Functions
changespecialchar(){
        # Function to change * characters in file lines from being read as wildcards when echo'd 
        # back outside of while read loop 
        
    line=$1
    # replace the wildcard symbol with a simple str of simply the * character
    line=${line//"*"/\"\*\"}
    # replacing quotation marks from multiple * so that full line is one single string
    searchfor1="\"\""
    searchfor2="\" "
    searchfor3=" \""
    line=${line//$searchfor1/}
    line=${line//$searchfor2/" "}
    line=${line//$searchfor3/" "}
}

obtainingerrorcode(){
        # Function to extract the error code of a given error line argument
        
    errorline=$1
    # searching for last instance of a double blank space and removing all characters before last 
    # instance from string 
    search1="  "
    rest1=${errorline##*$search1}
    # Finding index of last character before start of error code
    index1=$(( ${#errorline} - ${#rest1} - ${#search1} ))
    # finding the index of the first character after the end of the error code
    search2=" \*\""
    index2=$(( ${#errorline} - ${#search2} ))
    # producing a substring of only error code. 
    length=$(( $index2 - $index1 + 1 ))
    code=${errorline:$index1:$length}
    # echo $code
}

emailerrorlinesbody(){
    # check whether error code has lines before 1st stars of box to be include in 
    # email
    if [[ $(echo ${errorcodes[@]} | fgrep -w $2) ]];
    then
            # add 1st stars of box to body of email 
            body="${LINES[${abtls[$1]}-3]}\n" # ${LINES[${abtls[$j]}-4]}"
            # find and add lines of error code before 1st stars to body of email
            cnt=0
            while [[ ! ${LINES[${abtls[$1]}-3+$cnt]} == *"Possible matches"*  ]];
            do
                    ((cnt--))
                    body="${LINES[${abtls[$1]}-3+$cnt]}\n$body"
            done 
    else
    # add 1st stars of box to body of email 
    body="${LINES[${abtls[$1]}-3]}\n"
    fi 
    # add rest of error box to body of email. 
    for (( d=-2 ; d<7 ; d++)) ;
    do
            body="$body${LINES[${abtls[$1]}+d]}\n"
    done
    indx=${abtls[$1]}+6
    stars="${LINES[$indx]}"
    # continue to add lines to body of email until another set of starts, end of file or "possible
    # match for"
    c=1
    while [[ ! ${LINES[$indx+$c]} == *"Possible matches"*  ]] && [[ ! ${LINES[$indx+$c]} == $stars ]] && [[ $indx+$c -lt $numoflines ]];
    do
            body="$body${LINES[$indx+$c]}\n"
            ((c++))
    done
}

send_email_error(){
        # Function to send email to self about an error occuring. 
        
    echo " - sending email about error"
    # "project name run type error" to put in subject line of email.
    inst1=($(grep -i "Project name" $File))
    inst2=($(grep -i "Run type" $File))
    subject="${inst1[3]} ${inst2[3]} error"
    # give 1 instance of each different error within the block of latest errors stated 
    # in log file within body of email 
    after=$1
    # making arrays to keep all lines and the specific line numbers where "ABORT" is found
    declare -a LINES
    declare -a abtls
    i=-1
    while read -r line
    do 
            ((i++))
            # if line contains a '*' which maybe be interrupted as a wildcard
            if [[ "$line" == *"*"* ]];
            then
                    changespecialchar "$line"
            fi
            LINES[$i]=$line
            if [[ $i -gt $after ]] && [[ $line == *"ABORT"* ]];
            then
                    abtls+=(${i})
            fi
    done <$File  
    numoflines="${#LINES[@]}"
    numoferrors="${#abtls[@]}"
    declare -a found_error_codes
    e_body="Errors which have occurred in the last calculation run are\n\n "
    for (( j=0 ; j<$numoferrors ; j++ )) ; 
    do
            if [[ $j -eq 0 ]];
            then 
                    # extracting error code and add it to array for found error codes 
                    obtainingerrorcode "${LINES[${abtls[$j]}+5]}"
                    found_error_codes[0]=$code
                    emailerrorlinesbody $j $code
                    e_body="$e_body$body\n\n"
            else
                    obtainingerrorcode "${LINES[${abtls[$j]}+5]}"
                    if [[ ! $(echo ${found_error_codes[@]} | fgrep -w $code) ]];
                    then
                           found_error_codes+=(${code})
                           emailerrorlinesbody $j $code
                           e_body="$e_body$body\n\n"
                    fi
            fi
    done
    bash ~/send_message.sh "${subject}" "${e_body}"
}

send_email_finished(){
        # Function to send email to self about job finishing 
        
    echo " - sending email that job has finished to self"
    # "project name run type finished" to put in subject line of email.
    inst1=($(grep -i "Project name" $File))
    inst2=($(grep -i "Run type" $File))
    subject="${inst1[3]} ${inst2[3]} finished"
    echo $subject
        
    # give all warnings stated within log file within body of email. 
    body=($(grep -i "warn" $File))
    echo ${body[@]}
        
    bash ~/send_message.sh "${subject}" "${body[@]}"
}

prep_inp_4_restart(){
        # Function to prepare the inp file from being the first input sran file of the calculation 
        # to being an inp file which the calculation can be continued with, via the SCF being 
        # restarted, when it is used for the next srun.
        
    echo " - pareparing input file for next calc run"
    # check if the directory for the intermediate input/restart files doesn't already exists. 
    if [ ! -d $dir ];
    then
            # if it doesn't exist create the directory
            mkdir $dir
    fi
            # replace the string of $atomicSCF with the string of restartSCF within $inp file while 
            # producing a new file, ${inp}-e, which is the version of $inp before any strings were 
            # replaced
            sed -i -e "s/${atomic_SCF}/${restart_SCF}/g" $inp
            mv ${inp}-e $dir/beginning_$inp
}

change_scf(){
        # Function to change keyword within $restart for SCF GUESS from ATOMIC to RESTART
        
    echo " - changing SCF setting"
    # replace the string of $atomicSCF with the string of restartSCF within $restart file while 
    # producing a new file, ${restart}-e, which is the version of $restart before any strings were 
    # replaced
    sed -i -e "s/${atomic_SCF}/${restart_SCF}/g" $restart
    # move version of $restart before sed command implemented into directory for intermediate 
    # input/restart files
    mv ${restart}-e $dir/${restart}_1
}

rmline4restart(){
        # Function to remove keyword lines within the restart file which are not needed for next 
        # rerun of the calculation
        
    echo " - removing unneeded lines for next run"
    # set up array restart file name as multiple elements for multiple sed commands to be used on
    # the same file one after the other - sed commands could not be combined and the same named 
    # variable cannot be used for the file name of subsequent sed commands. 
    restart_=(${restart} ${restart} ${restart})
    # remove the lines containing $WFNname, $COORDformat, $COORDfile from $restart fil
    sed -i -e "/${WFNname}/d" ${restart_[0]}
    sed -i -e "/${COORDformat}/d" ${restart_[1]}
    sed -i -e "/${COORDfile}/d" ${restart_[2]}
}

prep_restrt_file(){
        # Function to prepare the first of the restart files for a restart and continuation of the
        # calculation from the standard CP2K autogenerated restart file, which contains the SCF 
        # setting of the previous inp file. 
    
    echo " - preparing restart file"
    change_scf
    rmline4restart
    # version of $restart file before last sed command was implemented is not needed therefore 
    # than be deleted. 
    rm ${restart}-e
}

saveintermediaterestarts(){
        # Function to give an appropriate name to the archive of the current intermediate restart 
        # file being made.
        
    echo " - giving an unique name to next before changed intermediate restart file to be archieved"
    counter=1
    saved=${restart}_${counter}
    until [ ! -f $dir/$saved ] ;
    do
            ((counter++))
            saved=${restart}_${counter}
   done
   mv ${restart}-e $dir/$saved
}

archivingxyzfile(){
        # Function to give an appropriate name to the archive of the current xyz file being made
        
    echo " - archiving current xyz file before calc run with new max. iteration settings"
    counter=1
    saved=${xyz_}_${counter}
    until [ ! -f $dir/$saved ] ;
    do
            ((counter++))
            saved=${xyz_}_${counter}
   done
   mv ${xyz_} $dir/$saved
}

correct_max_iter(){
        # Function to change key word value for GEO OPT MAX_ITER based upon whether the value of
        # the next rerun's calculation starting step is above a certain threshold ie 400.
        
    echo " - correct_max_iter in progress"
    max_iter=($(grep -i "MAX_ITER" $restart))
    step_strt_val=($(grep -i "STEP_START_VAL" $restart))
    #
    if [[ ${max_iter[1]} -eq 500 ]] && [[ ${step_strt_val[1]} -gt 400 ]];
    then
            replace="${max_iter[0]}  ${max_iter[1]}"
            with="${max_iter[0]}  800"
            echo $replace $with
            sed -i -e "s/${replace}/${with}/g" $restart
            saveintermediaterestarts
            archivingxyzfile
    fi
}

check_restrt_file_ready(){
        # Function to ensure that the current restart file being used as the input for the next  
        # rsun is properly set up for the next restart and continuation of the calculation. Checks 
        # whether the restart file has already been changed for use to continue the calculation 
        # from the standard CP2K autogenerated restart file or whether the restart file needs to be 
        # prepared still.
        
    echo " - checking restart file is ready for next calc submission"
    # if the restart file needs to be prepared still for next restart and calc continuation.
    if  grep -q $WFNname $restart && grep -q $COORDformat $restart && grep -q $COORDfile $restart && grep -q "$atomic_SCF" $restart;
    then 
            echo " - restart file still needs to be prepared"
            prep_restrt_file
    elif grep -q $WFNname $restart && grep -q $COORDformat $restart && grep -q $COORDfile $restart;
    then
            echo " - some lines need to be removed before next run"
            rmline4restart
            saveintermediaterestarts
    elif grep -q "$atomic_SCF" $restart;
    then 
            change_scf
    # if restart file has already been previously prepared.
    else
            echo " - restart file has been prepared previously, max iterations still to be checked"
            correct_max_iter
    fi
}

check_opt_step(){
        # Function to work out which optimization step the last run of the calculation ended on and 
        # direct the code execution to the appropriate stream of code based on the last  
        # optimization step. 
        
    echo " - checking optimization step"
    # save all instances found via grep of "OPTIMIZATION STEP:" within $File as an array
    instances=($(grep -i "OPTIMIZATION STEP:" $File))
    # find length of/number of elements within array
    length="${#instances[@]}"
    # store the last element of the list as this is the last optimization step to have run.
    step="${instances[$length-1]}"  
    echo "- last optimization step: $step"
    if [[ $step -eq 2 ]];
    then
            echo " - last optimization step run was step 2"
            prep_restrt_file
            ## srun --hint=nomultithread --distribution=block:block cp2k.psmp -i $restart >> $File &
    elif [[ $step -lt 2 ]];
    then
            echo " - picked up that less than 2 optimization steps have run for calculation"
            prep_inp_4_restart
            ## srun --hint=noultithread --distribution=block:block cp2k.psmp -i $inp >> $File &
    else
            echo " - restart file to be checked"
            check_restrt_file_ready
            ## srun --hint=nomultithread --distribution=block:block cp2k.psmp -i $restart >> $File &
    fi
}

check_4_ending_errors(){
        # Function to determine whether the job ending error identified by the grep command which 
        # directed the code execution to this function has occured in the last calculation run or 
        # previous to the last calculation run and has been fixed since occurring. 
        
    echo " - Investigating job ending error"
    i=0
    start=0
    error=0
    # 
    while read line
    do 
        ((i++))
        # if "ABORT" is found within the current $line at index $i, resign the value of error to be 
        # that of the line index.
        if [[ $line == *"ABORT"* ]];
        then
            error=$i 
        # if "PROGRAM STARTED AT" is found within the current $line at index $i, resign the value 
        # of start to be that of the line index.
        elif [[ $line == *"PROGRAM STARTED AT"* ]];
        then
            start=$i
        fi
    done <$File
    # if job ending error occurred within the last calculation run.
    if [[ $error -gt $start ]]; 
    then 
            echo " - Job ending error occured within last calculation run. Will send email to self about error"
            send_email_error $start
    # if error has been fixed and the calculation has been run since the job ending error occurred.
    else 
            echo " - error already fixed. Will want to go forward with other checks before srun command"
            check_opt_step
    fi
}

# check log file exists
if [[ -f $File ]];
then
        # if -c flag included in sbatch command, job ending errors do not need to be checked for. 
        if [[ $corrected == 'true' ]];
        then
                echo " - looking for job ending error codes skipped. Will want other checks to take place and srun command to be stated"
                check_opt_step
        # check log file for job ending error codes, job ending error found 
        elif grep -q "ABORT" $File ;
        then 
                echo " - job ending error has been found"
                check_4_ending_errors
        # 
        elif grep -q "PROGRAM ENDED AT" $File ;
        then
                echo " - Calculation has finished"
                # check if the directory for slurm files doesn't already exists. 
                if [ ! -d slurm_files ];
                then
                        # if it doesn't exist create the directory
                        mkdir slurm_files
                fi
                mv *.slurm slurm_files
                send_email_finished
        # no job ending errors found 
        else 
                echo " - no error found in file, other checks to take place"
                check_opt_step
        fi
else 
        # if log file doesn't exist then srun command should be run with .inp file of calculation
        echo " - $File does not exist yet, srun command to be run with .inp file, then resubmit_file function"
        ## srun --hint=noultithread --distribution=block:block cp2k.psmp -i $inp >> $File &
fi