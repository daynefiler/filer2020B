int[] tst = {0 , 0};
for (SAMRecord record : records) {
  if (record.getAttribute("XA") != null) tst[0]++;
  if (record.getMappingQuality() < 30) tst[1]++;
}
return tst[0] < 2 && tst[1] == 0; 
