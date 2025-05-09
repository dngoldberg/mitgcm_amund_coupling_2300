for i in run_ice_2009_TC_80_iceParm3_coul; do 
	cd $i; 
	for j in meta data; do 
		for k in pickup pickup_streamice; do 
			ls -lh $k.ckptA_old.$j 
			cp $k.ckptA_old.$j $k.ckptA.$j -v; 
		done;
#	        cp pickup_streamice.ckptA.old.$j pickup_streamice.ckptAlast.$j -v	
		for k in pickup_streamice; do 
			ls -lh $k.ckptA_old.$j
			cp $k.ckptA_old.$j $k.ckptAlast.$j -v; 
		done;
	done; 
	cd ..; 
done
