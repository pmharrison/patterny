# 
# AWK program, part of patterny shell script 
# 
# -v prefix 
# 
# 
BEGIN{flag=0;} 
NR>0{flag=1;}
END{if(flag){print "Patterny warning: the file prefix " prefix " is already populated in this directory...Please make sure it is specific enough...\n";}} 

