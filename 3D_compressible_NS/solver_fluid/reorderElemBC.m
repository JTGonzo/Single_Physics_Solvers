function [cnnBCnew] = reorderElemBC(BC, cnnBC, Flag3D)
	if(~Flag3D)
		for i=1:size(BC,1)
		   [tf, idx1(i)] = ismember(BC(i,1),cnnBC(i,:));
		   cnnBCnew(i,1) = cnnBC(i,idx1(i));
		   if (idx1(i) == 1)
			   cnnBCnew(i,2) = cnnBC(i,2);
			   cnnBCnew(i,3) = cnnBC(i,3);
			   cnnBCnew(i,4) = cnnBC(i,4);    
		   elseif (idx1(i) == 2)
			   cnnBCnew(i,2) = cnnBC(i,3);
			   cnnBCnew(i,3) = cnnBC(i,4);
			   cnnBCnew(i,4) = cnnBC(i,1);
		   elseif (idx1(i) == 3)
			   cnnBCnew(i,2) = cnnBC(i,4);
			   cnnBCnew(i,3) = cnnBC(i,1);
			   cnnBCnew(i,4) = cnnBC(i,2);        
		   elseif (idx1(i) == 4)
			   cnnBCnew(i,2) = cnnBC(i,1);
			   cnnBCnew(i,3) = cnnBC(i,2);
			   cnnBCnew(i,4) = cnnBC(i,3);
		   end
		   [tf, idx2(i)] = ismember(BC(i,2),cnnBCnew(i,:));
		   if (idx2(i) == 4)
			   temp = cnnBCnew(i,1);
			   cnnBCnew(i,1) = cnnBCnew(i,4);
			   cnnBCnew(i,4) = cnnBCnew(i,3);
			   cnnBCnew(i,3) = cnnBCnew(i,2);
			   cnnBCnew(i,2) = temp;
           end		   
		   clear tf idx1 idx2		   
		end
	else
		for i=1:size(BC,1)
		   [tf, idx1(i)] = ismember(BC(i,1),cnnBC(i,:));
		   cnnBCnew(i,1) = cnnBC(i,idx1(i));
		   if (idx1(i) == 1)
			   cnnBCnew(i,2) = cnnBC(i,2);
			   cnnBCnew(i,3) = cnnBC(i,3);
			   cnnBCnew(i,4) = cnnBC(i,4); 
			   cnnBCnew(i,5) = cnnBC(i,5);
			   cnnBCnew(i,6) = cnnBC(i,6);
			   cnnBCnew(i,7) = cnnBC(i,7);
			   cnnBCnew(i,8) = cnnBC(i,8);			   
		   elseif (idx1(i) == 2)
			   cnnBCnew(i,2) = cnnBC(i,3);
			   cnnBCnew(i,3) = cnnBC(i,4);
			   cnnBCnew(i,4) = cnnBC(i,5);
			   cnnBCnew(i,5) = cnnBC(i,6);
			   cnnBCnew(i,6) = cnnBC(i,7);
			   cnnBCnew(i,7) = cnnBC(i,8);
			   cnnBCnew(i,8) = cnnBC(i,1);
		   elseif (idx1(i) == 3)
			   cnnBCnew(i,2) = cnnBC(i,4);
			   cnnBCnew(i,3) = cnnBC(i,5);
			   cnnBCnew(i,4) = cnnBC(i,6);
			   cnnBCnew(i,5) = cnnBC(i,7);
			   cnnBCnew(i,6) = cnnBC(i,8);
			   cnnBCnew(i,7) = cnnBC(i,1);
			   cnnBCnew(i,8) = cnnBC(i,2);			   
		   elseif (idx1(i) == 4)
			   cnnBCnew(i,2) = cnnBC(i,5);
			   cnnBCnew(i,3) = cnnBC(i,6);
			   cnnBCnew(i,4) = cnnBC(i,7);
			   cnnBCnew(i,5) = cnnBC(i,8);
			   cnnBCnew(i,6) = cnnBC(i,1);
			   cnnBCnew(i,7) = cnnBC(i,2);
			   cnnBCnew(i,8) = cnnBC(i,3);		
		   elseif (idx1(i) == 5)
			   cnnBCnew(i,2) = cnnBC(i,6);
			   cnnBCnew(i,3) = cnnBC(i,7);
			   cnnBCnew(i,4) = cnnBC(i,8);
			   cnnBCnew(i,5) = cnnBC(i,1);
			   cnnBCnew(i,6) = cnnBC(i,2);
			   cnnBCnew(i,7) = cnnBC(i,3);
			   cnnBCnew(i,8) = cnnBC(i,4);
		   elseif (idx1(i) == 6)
			   cnnBCnew(i,2) = cnnBC(i,7);
			   cnnBCnew(i,3) = cnnBC(i,8);
			   cnnBCnew(i,4) = cnnBC(i,1);
			   cnnBCnew(i,5) = cnnBC(i,2);
			   cnnBCnew(i,6) = cnnBC(i,3);
			   cnnBCnew(i,7) = cnnBC(i,4);
			   cnnBCnew(i,8) = cnnBC(i,5);
		   elseif (idx1(i) == 7)
			   cnnBCnew(i,2) = cnnBC(i,8);
			   cnnBCnew(i,3) = cnnBC(i,1);
			   cnnBCnew(i,4) = cnnBC(i,2);
			   cnnBCnew(i,5) = cnnBC(i,3);
			   cnnBCnew(i,6) = cnnBC(i,4);
			   cnnBCnew(i,7) = cnnBC(i,5);
			   cnnBCnew(i,8) = cnnBC(i,6);	
		   elseif (idx1(i) == 8)
			   cnnBCnew(i,2) = cnnBC(i,1);
			   cnnBCnew(i,3) = cnnBC(i,2);
			   cnnBCnew(i,4) = cnnBC(i,3);
			   cnnBCnew(i,5) = cnnBC(i,4);
			   cnnBCnew(i,6) = cnnBC(i,5);
			   cnnBCnew(i,7) = cnnBC(i,6);
			   cnnBCnew(i,8) = cnnBC(i,7);			   
		   end
		   [tf, idx2(i)] = ismember(BC(i,2),cnnBCnew(i,:));
		   if (idx2(i) == 8)
			   temp = cnnBCnew(i,1);
			   cnnBCnew(i,1) = cnnBCnew(i,8);
			   cnnBCnew(i,8) = cnnBCnew(i,7);
			   cnnBCnew(i,7) = cnnBCnew(i,6);
			   cnnBCnew(i,6) = cnnBCnew(i,5);
			   cnnBCnew(i,5) = cnnBCnew(i,4);
			   cnnBCnew(i,4) = cnnBCnew(i,3);		   
			   cnnBCnew(i,3) = cnnBCnew(i,2);
			   cnnBCnew(i,2) = temp;
		   end
		   [tf, idx3(i)] = ismember(BC(i,3),cnnBCnew(i,:));
		   if (idx3(i) == 8)
			   temp = cnnBCnew(i,1);
			   cnnBCnew(i,1) = cnnBCnew(i,8);
			   cnnBCnew(i,8) = cnnBCnew(i,7);
			   cnnBCnew(i,7) = cnnBCnew(i,6);
			   cnnBCnew(i,6) = cnnBCnew(i,5);
			   cnnBCnew(i,5) = cnnBCnew(i,4);
			   cnnBCnew(i,4) = cnnBCnew(i,3);		   
			   cnnBCnew(i,3) = cnnBCnew(i,2);
			   cnnBCnew(i,2) = temp;
		   end		
		   [tf, idx4(i)] = ismember(BC(i,4),cnnBCnew(i,:));
		   if (idx4(i) == 8)
			   temp = cnnBCnew(i,1);
			   cnnBCnew(i,1) = cnnBCnew(i,8);
			   cnnBCnew(i,8) = cnnBCnew(i,7);
			   cnnBCnew(i,7) = cnnBCnew(i,6);
			   cnnBCnew(i,6) = cnnBCnew(i,5);
			   cnnBCnew(i,5) = cnnBCnew(i,4);
			   cnnBCnew(i,4) = cnnBCnew(i,3);		   
			   cnnBCnew(i,3) = cnnBCnew(i,2);
			   cnnBCnew(i,2) = temp;
		   end		
		   clear tf idx1 idx2 idx3 idx4
		end
    end
end