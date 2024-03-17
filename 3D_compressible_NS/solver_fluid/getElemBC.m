function [cnnBC] = getElemBC(BC, cnn, Flag3D)

	[tf,idx]      = ismember(BC(:,1),cnn(:,1));
	[tf,idx(:,2)] = ismember(BC(:,1),cnn(:,2));
	[tf,idx(:,3)] = ismember(BC(:,1),cnn(:,3));
	[tf,idx(:,4)] = ismember(BC(:,1),cnn(:,4));
    
	if(Flag3D)
		[tf,idx(:,5)] = ismember(BC(:,1),cnn(:,5));
		[tf,idx(:,6)] = ismember(BC(:,1),cnn(:,6));
		[tf,idx(:,7)] = ismember(BC(:,1),cnn(:,7));	
		[tf,idx(:,8)] = ismember(BC(:,1),cnn(:,8));	
    end   
	elem = unique(idx(:), 'stable');
% 	elem = unique(idx(:));
	elem(elem==0)=[];

	[tf,idx2]      = ismember(BC(:,2),cnn(:,1));
	[tf,idx2(:,2)] = ismember(BC(:,2),cnn(:,2));
	[tf,idx2(:,3)] = ismember(BC(:,2),cnn(:,3));
	[tf,idx2(:,4)] = ismember(BC(:,2),cnn(:,4));
    
	if(Flag3D)
		[tf,idx2(:,5)] = ismember(BC(:,2),cnn(:,5));
		[tf,idx2(:,6)] = ismember(BC(:,2),cnn(:,6));
		[tf,idx2(:,7)] = ismember(BC(:,2),cnn(:,7));	
		[tf,idx2(:,8)] = ismember(BC(:,2),cnn(:,8));		
    end	    
	elem2 = unique(idx2(:), 'stable');
% 	elem2 = unique(idx2(:));
	elem2(elem2==0)=[];
	
	if(Flag3D)
		[tf,idx3]      = ismember(BC(:,3),cnn(:,1));
		[tf,idx3(:,2)] = ismember(BC(:,3),cnn(:,2));
		[tf,idx3(:,3)] = ismember(BC(:,3),cnn(:,3));
		[tf,idx3(:,4)] = ismember(BC(:,3),cnn(:,4));
		[tf,idx3(:,5)] = ismember(BC(:,3),cnn(:,5));
		[tf,idx3(:,6)] = ismember(BC(:,3),cnn(:,6));
		[tf,idx3(:,7)] = ismember(BC(:,3),cnn(:,7));	
		[tf,idx3(:,8)] = ismember(BC(:,3),cnn(:,8));		
		elem3 = unique(idx3(:), 'stable');
% 		elem3 = unique(idx3(:));
		elem3(elem3==0)=[];
		
		[tf,idx4]      = ismember(BC(:,4),cnn(:,1));
		[tf,idx4(:,2)] = ismember(BC(:,4),cnn(:,2));
		[tf,idx4(:,3)] = ismember(BC(:,4),cnn(:,3));
		[tf,idx4(:,4)] = ismember(BC(:,4),cnn(:,4));
		[tf,idx4(:,5)] = ismember(BC(:,4),cnn(:,5));
		[tf,idx4(:,6)] = ismember(BC(:,4),cnn(:,6));
		[tf,idx4(:,7)] = ismember(BC(:,4),cnn(:,7));	
		[tf,idx4(:,8)] = ismember(BC(:,4),cnn(:,8));		
		elem4 = unique(idx4(:), 'stable');
% 		elem4 = unique(idx4(:));
		elem4(elem4==0)=[];				
    end    
    elem = intersect(elem, elem2, 'stable');

	if(Flag3D)
		elem = intersect(elem, elem3, 'stable');
        elem = intersect(elem, elem4, 'stable');
	end

	cnnBC = cnn(elem,:);
    
	clear elem elem2 idx idx2
	if(Flag3D)
		clear elem3 elem4 idx3 idx4
    end
end