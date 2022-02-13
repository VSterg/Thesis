function [C CM CM_metrics D] = compute_network_variance(A,B)
%The purpose of this function is to compute the difference between two
%networks by comparing their mixing matrices.
%Inputs: 
%A : network (ground truth if available) mixing matrix
%B : network (estimated) mixing matrix 
%Outputs:
%C : number of missclassified connections (not present in ground truth
%network)
%D : error matrix
%CM : confusion matrix, CM(1,1): correctly identified as not connected
%                       CM(1,2): incorrectly identified as connected
%                       CM(2,1): incorrectly identified as not connected
%                       CM(2,2): correctly identified as connected
%missclassification_matrix : logical matrix where 1 = missclassification
C = A - B;
C = C.*C;
missclassification_matrix = C;
D = C;
C = sum(sum(C)); 
A = reshape(A,[1 size(A,1)^2]);
B = reshape(B,[1 size(B,1)^2]);
A = logical(A);
B = logical(B);
CM = confusionmat(A,B); 
CM_metrics = compute_conf_mat_metrics(CM);
