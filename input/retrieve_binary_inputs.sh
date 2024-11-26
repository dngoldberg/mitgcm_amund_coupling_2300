input_dir=/home/dgoldber/network_links/geosIceOcean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling/input

for i in 2009;
do
 echo "getting all ${i}"
 for j in ice oce;
 do
  echo start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/*.bin ./start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/*.box ./start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/*obw* ./start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/*obs* ./start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/*.init.* ./start_${i}_input/$j
  scp -P 6022 -i /home/n02/n02/dngoldbe/.ssh/id_rsa_cirrus_geos dgoldber@sshpa.geos.ed.ac.uk:${input_dir}/start_${i}_input/$j/pickup* ./start_${i}_input/$j
 done
done

