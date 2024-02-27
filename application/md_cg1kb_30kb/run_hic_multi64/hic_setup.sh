awk '
 BEGIN{
  Nbp=1000 # bps per bead,size of the bin
  N0=31177000
  bins=30
  for(i=1;i<=bins;i+=1){
   for(j=i;j<=bins;j+=1){
    reads[i,j]=0.0
   }
  }
 }
 {
  i=1+($1-N0)/Nbp
  j=1+($2-N0)/Nbp
  reads[i,j]=$3
 }
 END{
  for(i=1;i<=bins;i+=1){
   for(j=i;j<=bins;j+=1){
    print i,j,reads[i,j]
   }
  }
 }
' ../../experimental_data/hic_chr5-31177kb-31207kb_1000bp_coverage.dat > hic_data_all.dat
r0=1
d0=0.5
awk '
 BEGIN{
  print "HIC2 ..."
  print "R_0="0.1*'$r0' # the unit is nm
  print "D_0="0.1*'$d0' # in nm
  i=0
  delta_min=1 # minimum spacing to include a Hi-C contact
 }
 {
  p1=$1
  p2=$2
  reads=$3
  delta=p2-p1
  strength_min=0.0
  if(delta>=delta_min){
   i=i+1
   print "ATOMS"i"="p1","p2" REFERENCE"i"="reads" # "p1,p2
  }
 }
 END{
  print "LABEL=hic"
  print "... HIC2"
 }
' hic_data_all.dat > plumed_hic_all.dat
