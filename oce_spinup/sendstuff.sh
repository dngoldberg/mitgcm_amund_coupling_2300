for yr in 1995 2009; do
 for suff in data meta; do
   for file in pickup pickup_streamice pickup_shelfice; do

  scp -i ~/.ssh/id_rsa_cirrus_geos run_${yr}_1_851/$file.ckptB.$suff dgoldber@ssh.geos.ed.ac.uk:/home/dgoldber/network_links/iceOceanShare/dgoldber/archer_output/AMUND_COUPLE/OCE_SPINUP/run_${yr}_1_851/$file.0000933120.$suff

 done
done
done
