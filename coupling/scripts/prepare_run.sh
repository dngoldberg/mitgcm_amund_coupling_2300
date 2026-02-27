if [ -d "../run_ice_$1_$2_$3_$4_$5" ]; then
  cd ../run_ice_$1_$2_$3_$4_$5
  rm -rfv *
  mkdir diags
else
  echo 'There is no ice run directory'
  mkdir ../run_ice_$1_$2_$3_$4_$5
  mkdir ../run_ice_$1_$2_$3_$4_$5/diags
fi

if [ -d "../run_oce_$1_$2_$3_$4_$5" ]; then
  cd ../run_oce_$1_$2_$3_$4_$5
  rm -rfv *
  mkdir diags
else
  echo 'There is no oce run directory'
  mkdir ../run_oce_$1_$2_$3_$4_$5
  cd ../run_oce_$1_$2_$3_$4_$5
  mkdir diags
fi
cp -f ../../input/start_$1_input/oce/pickup* ./


if [ -d "../stdout_oce_files_$1_$2_$3_$4_$5" ]; then
  cd ../stdout_oce_files_$1_$2_$3_$4_$5
  rm -rfv *
else
  echo 'There is no ice run directory'
  mkdir ../stdout_oce_files_$1_$2_$3_$4_$5
fi

if [ -d "../stdout_ice_files_$1_$2_$3_$4_$5" ]; then
  cd ../stdout_ice_files_$1_$2_$3_$4_$5
  rm -rfv *
else
  echo 'There is no ice run directory'
  mkdir ../stdout_ice_files_$1_$2_$3_$4_$5
fi

if [ -d "../transfer_files_$1_$2_$3_$4_$5" ]; then
  cd ../transfer_files_$1_$2_$3_$4_$5
  rm -rfv *
else
  echo 'There is no ice run directory'
  mkdir ../transfer_files_$1_$2_$3_$4_$5
fi

echo "DONE WITH PREP"
