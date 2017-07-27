function th=arc_angles(ang,d)
if d>1
    d=1;
end
if d<-1
    d=-1;
end

if strcmp('acos',ang)==1
    th=acos(d);
end
if strcmp('asin',ang)==1
    th=asin(d);
end  

if th<0
    th=th+2*pi;
end
if th>2*pi
    th=th-2*pi;
end
    

