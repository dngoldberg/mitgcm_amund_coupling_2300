#!/bin/bash
################################################
# Clean out old results and link input files.
################################################

rundir=run_$1_$2_$3_$4_$5_$6

# $1 year
# $2 PAS
# $3 G/NG
# $4 cfric
# $5 heattranscoeff
# $6 dig

echo "rundir"
echo $rundir
# Empty the run directory - but first make sure it exists!
if [ -d "../$rundir" ]; then
  cd ../$rundir
  rm -rf *
else
  echo 'creating run directory'
  mkdir ../$rundir
  cd ../$rundir
fi

# Link everything from the input directory
input_oce_dir=../../input/start_$1_input/oce
bmyear=''

cd $input_oce_dir
bash extend_obcs_back.sh 
cd $OLDPWD

ln -s $input_oce_dir/* . 


# Deep copy of the master namelist (so it doesn't get overwritten in input/)
rm -f data
cp -f $input_oce_dir/data .
rm -f data.shelfice
cp -f $input_oce_dir/data.shelfice .

if [[ $1 == 2004 ]]; then 
 yrfile=2001
 calstring=" startdate_1 = 20010101,"
else 
 yrfile=2008
 calstring=" startdate_1 = 20080101,"
fi

rm data.cal
cp -f $input_oce_dir/data.cal .
sed "s/.*startdate_1.*/$calstring/" data.cal > data.temp;
mv data.temp data.cal

newdata=" shiCdrag = $4,"
sed "s/.*shiCdrag.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newobcs=" obsufile = 'uvel.obs.$2'",
sed "s/.*obsufile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obsvfile = 'vvel.obs.$2'",
sed "s/.*obsvfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obssfile = 'salt.obs.$2'",
sed "s/.*obssfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obstfile = 'temp.obs.$2'",
sed "s/.*obstfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwufile = 'uvel.obw.$2'",
sed "s/.*obwufile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwvfile = 'vvel.obw.$2'",
sed "s/.*obwvfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwsfile = 'salt.obw.$2'",
sed "s/.*obwsfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs
newobcs=" obwtfile = 'temp.obw.$2'",
sed "s/.*obwtfile.*/$newobcs/" data.obcs > data.temp; mv data.temp data.obcs

if [ $3 == G ]; then
	newdata=" SHELFICEuseGammaFrict = .true.,"
else
	newdata=" SHELFICEuseGammaFrict = .false.,"
fi
sed "s/.*useGammaFrict.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice

newdata=" psurfinitfile = 'etainit.$2.bin',"
sed "s/.*psurfinitfile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" shelficetopofile = 'shelftopo.$2.bin',"
sed "s/.*shelficetopofile.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newdata=" hydrogsaltfile = 'salt.init.$yrfile.$2'"
sed "s/.*hydrogsaltfile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" hydrogthetafile = 'theta.init.$yrfile.$2'"
sed "s/.*hydrogthetafile.*/$newdata/" data > data.temp; mv data.temp data

if [ $# == 6 ]; then
newdata=" bathyfile = 'bathy_mod_dig_$6.bin',"
else
newdata=" bathyfile = 'bathy_mod.bin',"
fi
sed "s/.*bathyfile.*/$newdata/" data > data.temp; mv data.temp data
newdata=" shelficemassfile = 'shelficemassinit.bin',"
sed "s/.*shelficemassfile.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
newdata=" SHELFICEMassStepping = .true.,"
sed "s/.*SHELFICEMassStepping.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice

if [ $# == 6 ]; then
newdata=" streamicetopogfile = 'bathy_init_oce_dig_$6.bin',"
else	
newdata=" streamicetopogfile = 'bathy_init_oce.bin',"
fi
sed "s/.*streamicetopogfile.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice
newdata=" streamicethickfile = 'icethick_oce.bin',"
sed "s/.*streamicethickfile.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice
newdata=" streamicehmaskfile = 'hmask_oce.bin',"
sed "s/.*streamicehmaskfile.*/$newdata/" data.streamice > data.temp; mv data.temp data.streamice

#SHELFICE_transition_gamma
#SHELFICETransGammaThickness
# SHELFICEheatTransCoeff

if [ $3 != G ] && [ $3 != NG ]; then
	newdata=" SHELFICE_transition_gamma=.true.,"
	sed "s/.*SHELFICE_transition_gamma.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
	newdata=" SHELFICETransGammaThickness=$3,"
	sed "s/.*SHELFICETransGammaThickness.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
	newdata=" SHELFICEuseGammaFrict = .true.,"
	sed "s/.*useGammaFrict.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
fi

if [ $3 != G ]; then
	newdata=" SHELFICEheatTransCoeff=$5,"
        sed "s/.*SHELFICEheatTransCoeff.*/$newdata/" data.shelfice > data.temp; mv data.temp data.shelfice
fi



# Deep copy of any pickups (so they don't get overwritten in input/)
rm -f pickup*
#cp -f $input_oce_dir/pickup* . 2>/dev/null

# Link forcing files stored elsewhere

# Link executables
ln -s ../build/mitgcmuv .
ln -s ../../python_scripts/* .
#ln -s ../../../utilities/mit2nc/mit2nc .
