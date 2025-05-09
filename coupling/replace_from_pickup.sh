init=TC
icedir=run_ice_2009_${init}_80_iceParm3_coul
ocedir=run_oce_2009_${init}_80_iceParm3_coul

for dir in $icedir $ocedir; do

for ext in data meta; do

if [ -z $1 ]; then
 latest_pickup_streamice=$(ls -t ${dir}/pickups/pickup_streamice*.${ext} 2>/dev/null | head -n 1)
 latest_pickup=$(ls -t ${dir}/pickups/pickup.*.${ext} 2>/dev/null | head -n 1)
else
 ocenum=$1
 icenum=$(( ($ocenum-1555200) / 8640 ))

 if [ $dir == $ocedir ]; then
  num=$ocenum
 else
  num=$icenum
 fi

 latest_pickup_streamice=$(ls -t ${dir}/pickups/pickup_streamice*${num}.${ext} 2>/dev/null | head -n 1)
 latest_pickup=$(ls -t ${dir}/pickups/pickup.*${num}.${ext} 2>/dev/null | head -n 1)
fi



cp $latest_pickup_streamice $dir/pickup_streamice.ckptA.${ext} -v
cp $latest_pickup $dir/pickup.ckptA.${ext} -v

if [ $dir == $ocedir ]; then
	if [ -z "$1" ]; then
	  num=$ocenum
  	  latest_pickup_shelfice=$(ls -t ${dir}/pickups/pickup_shelfice.${ext} 2>/dev/null | head -n 1)
  	else
  	  latest_pickup_shelfice=$(ls -t ${dir}/pickups/pickup_shelfice*${num}.${ext} 2>/dev/null | head -n 1)
	fi
	cp $latest_pickup_shelfice $dir/pickup_shelfice.ckptA.${ext} -v
	
fi

if [ $dir == $icedir ]; then
	if [ -z "$1" ]; then 
	  num=$icenum
	  latest_pickup_streamice=$(ls -t ${dir}/pickups/pickup_streamice*.${ext} 2>/dev/null | head -n 1)
	else
          latest_pickup_streamice=$(ls -t ${dir}/pickups/pickup_streamice*${num}.${ext} 2>/dev/null | head -n 1)
	fi
	cp $latest_pickup_streamice $dir/pickup_streamice.ckptAlast.${ext} -v
fi

done
done

#for i in run_ice_2009_TC_80_iceParm3_coul; do 
#	cd $i; 
#	for j in meta data; do 
#		for k in pickup pickup_streamice; do 
#			ls -lh $k.ckptA_old.$j 
#			cp $k.ckptA_old.$j $k.ckptA.$j -v; 
#		done;
#	        cp pickup_streamice.ckptA.old.$j pickup_streamice.ckptAlast.$j -v	
#		for k in pickup_streamice; do 
#			ls -lh $k.ckptA_old.$j
#			cp $k.ckptA_old.$j $k.ckptAlast.$j -v; 
#		done;
#	done; 
#	cd ..; 
#done
