function CSM_export_solution(dim, u, vertices, elements, nov, outputFileName, iter, variableName)

titleData = 'CSM_solution';
if nargin < 7
      iter = -1;
end

if nargin < 8
      variableName = 'StructureDisplacement';
end

novP1 = size(vertices,2);
      
      data=struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:4,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'VECTORS'}},...
            'variableData',{{[full(u(1:novP1)) full(u(nov+1:novP1+nov))  full(u(2*nov+1:novP1+2*nov))]'}});
         
timewrite = tic;
%
%  Get the size of your mesh
%
[ nodeNum, dim ] = size ( data.vertices );
[ elementNum,  nodePerElement ] = size ( data.elements );
%
%  Open the output file.
%
if data.iteration == -1
    output_filename = sprintf('%s.vtk', data.outputFile);
else
    output_filename = sprintf('%s%04d.vtk', data.outputFile,data.iteration);
end

if ( isempty ( output_filename ) )
    output_filename = 'solution.vtk';
end

outputSolution = fopen ( output_filename, 'w' );
%
%  Transpose or otherwise modify the data.
%
coorz = zeros ( 3, nodeNum );
coorz(1:dim,1:nodeNum) = (data.vertices(1:nodeNum,1:dim))';


elementNode = zeros ( nodePerElement, elementNum );
elementNode(1:nodePerElement,1:elementNum) = ...
    data.elements(1:elementNum,1:nodePerElement)' - 1;

vtk_headerfile ( outputSolution, data.title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode);
for i=1:size(data.variableName,2)
    vtk_write ( outputSolution, nodeNum, data.variableName{i}, data.variableType{i}, data.variableData{i});
end

fclose ( outputSolution );

timewrite = toc(timewrite);
fprintf ( 1, '\n' );
fprintf ( 1, '  The data was written to "%s" in %1.3f seconds\n', output_filename, timewrite);

return
end

function vtk_headerfile ( outputSolution, title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode)
%*****************************************************************************80
%
%% VTK_HEADERFILE write the header of the VTK file.

fprintf ( outputSolution, '# vtk DataFile Version 3.0\n' );
fprintf ( outputSolution, '%s\n', title );
fprintf ( outputSolution, 'BINARY\n\n' );
fprintf ( outputSolution, 'DATASET UNSTRUCTURED_GRID\n' );

fprintf ( outputSolution, 'POINTS %d float\n', nodeNum );
fwrite(outputSolution, coorz(1:3,:),'float','b');

cell_size = elementNum * ( nodePerElement + 1 );

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, '\nCELLS  %d  %d\n', elementNum, cell_size );


fwrite(outputSolution, [nodePerElement*ones(1,elementNum); elementNode],'int','b');


fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, '\nCELL_TYPES %d\n', elementNum );
fwrite(outputSolution, [10*ones(elementNum,1)]','int','b');


fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, '\nPOINT_DATA %d\n', nodeNum );
return
end

function vtk_write ( outputSolution, nodeNum, variableName, variableType, variableData)

if(variableType=='SCALARS')
    fprintf ( outputSolution, 'SCALARS %s float\n', variableName );
    fprintf ( outputSolution, 'LOOKUP_TABLE default\n' );
    
    fwrite(outputSolution, variableData(1,:),'float','b');
    fprintf ( outputSolution, '\n' );
    
else
    fprintf ( outputSolution, 'VECTORS %s float\n',variableName );
    
    fwrite(outputSolution, variableData(1:3,:),'float','b');
    fprintf ( outputSolution, '\n' );
end
return
end

