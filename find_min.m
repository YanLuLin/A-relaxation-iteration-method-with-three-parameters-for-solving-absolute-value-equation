function [min_element,location]=find_min(A)
m=min(A);          
min_element=min(m);          
[row,column]=find(A==min_element);         
location=[row,column];

