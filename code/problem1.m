%Matlab codes for Finite Element Analysis
clc; clear all;

elementNodes = [1 2;2 3;2 4];
numberElements = size(elementNodes,1);
numberNodes = 4;
displacements = zeros(numberNodes,1);
force = zeros(numberNodes,1);
stiffness = zeros(numberNodes);

force(2) = 10.0;
 
for e = 1:numberElements
    elementDof = elementNodes(e,:);
    stiffness(elementDof,elementDof) = ...
        stiffness(elementDof,elementDof) + [1 -1;-1 1];
end

%boundary conditions and solution
prescribedDof = [1;3;4];

%free Dof: activeDof
activeDof = setdiff([1:numberNodes]',[prescribedDof]);

%solution
%K U = F; U = K\F or U = F*inv(K)
displacements_tmp = stiffness(activeDof,activeDof)\force(activeDof);

displacements(activeDof)=displacements_tmp;
force = stiffness*displacements;

