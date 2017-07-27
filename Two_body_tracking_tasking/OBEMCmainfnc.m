function OBEMCmainfnc()

MOM1=cell(50,1);
MOM2=cell(50,1);
MOM3=cell(50,1);
MOM4=cell(50,1);

parfor i=1:50
    [M1,M2,M3,M4]=ORBmomsMC(i);
    MOM1{i}=M1;
    MOM2{i}=M2;
    MOM3{i}=M3;
    MOM4{i}=M4;
    
end
M1=0;
M2=0;
M3=0;
M4=0;
for i=1:1:50
  M1=M1+MOM1{i};
  M2=M2+MOM2{i};
  M3=M3+MOM3{i};
  M4=M4+MOM4{i};
    
end
save('OBEMCmoms50x100000','M1','M2','M3','M4')


end
