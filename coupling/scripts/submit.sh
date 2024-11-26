#!/bin/bash
################################################
# Start a self-resubmitting simulation.
################################################

# ID number for run

# clean run directory and link all required files

# record start times
TIMEQSTART="$(date +%s)"
echo Start-time `date` >> ../run_forward/times


# Initialize variables with default values
r_flag=false

# Function to display usage information
usage() {
    echo "Usage: $0 -y <value> -c <value> -a <value> -p <value> [-n <parm_value>] [-f <cfric_value>] [-s <shitrans_value>] [-d <depthtrans_value>] [-i <iceiter>] [-r]" 
    exit 1
}

# submitted Submitted batch job 7075393 -y 2009 -c TC -a 20 -p n02-GRISLAKES -n iceParm4 -f 0.014 -s 0.00014 -d G -i 25 -r false


r_flag=false
f_flag=false
s_flag=false
n_flag=false
d_flag=false
i_flag=false

# Parse command line options
while getopts ":y:c:a:p:f:s:r" opt; do
    case $opt in
        y)
            y_value="$OPTARG"
            ;;
        c)
            cal_value="$OPTARG"
            ;;
        a)
            PAS_value="$OPTARG"
            ;;
        p)
            p_value="$OPTARG"
            ;;
        r)
            r_flag=true
            ;;
	f)
            cfric="$OPTARG"
	    f_flag=true
	    ;;
	s)
            shitrans="$OPTARG"
	    s_flag=true
	    ;;
        n)
            parmn="$OPTARG"
	    n_flag=true
	    ;;
        d)
	    depth="$OPTARG"
	    d_flag=true
	    ;;
        i)
	    iter="$OPTARG"
            i_flag=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done



if [ $s_flag == true ]; then
  echo "using alternate shitrans"
else
  echo "using default shitrans"
  shitrans=0.00014
fi

if [ $f_flag == true ]; then
  echo "using alternate cfric"
else
  echo "using default cfric"
  cfric=0.01
fi

if [ $n_flag == true ]; then
  echo "using alternate parmname"
else
  echo "using default parmname"
  parmn='iceParmDefault'
fi

if [ $d_flag == true ]; then
  echo "using alternate depth code"
else
  echo "using default depth code"
  depth='G'
fi

if [ $i_flag == true ]; then
  echo "using alternate iter"
else
  echo "using default iter"
  iter=0
fi

if [ $r_flag == true ]; then
  echo "restarting."
else
  echo "The -r flag is not set. deleting files"
  bash prepare_run.sh $y_value $cal_value $PAS_value $parmn coul
fi

add_leading_zero() {
    # Check if the number starts with a dot
    if [[ $1 == .* ]]; then
        echo "0$1"
    else
        echo "$1"
    fi
}

cfric=$(add_leading_zero $cfric)
shitrans=$(add_leading_zero $shitrans)

#$y_value
#$cal_value
#$PAS_value
# coul
#20
#$p_value
#$cfric
#$shitrans
#$shitransdepth
#$n_value
#$iterTC

echo $TIMEQSTART
echo $p_value
# submit the job chain
echo "sbatch --job-name=c${cal_value}${PAS_value}$parmn -A $p_value run_repeat.slurm $y_value $cal_value $PAS_value coul 20 $p_value $cfric $shitrans $depth $parmn $iter" 
RES=$(sbatch --job-name=c${cal_value}${PAS_value}$parmn -A $p_value run_repeat.slurm $y_value $cal_value $PAS_value coul 20 $p_value $cfric $shitrans $depth $parmn $iter)

echo submitted $RES -y ${y_value} -c ${cal_value} -a ${PAS_value} -p ${p_value} -n $parmn -f $cfric -s $shitrans -d $depth -i $iter -r $r_flag >> job_id_list
