function [Vsort,Dsort, fp_modes] = getEigsAndSortModes(T0)
    global Ny;
    
    [Vector, sort_me] = eig(T0);
    
%sort_me = diag(D); %turn diagonal into column

s_fp = [0;1];
s_fe = [0;1];
s_bp = [0;1];
s_be = [0;1];

fp_j = 1;
fe_j = 1;
bp_j = 1;
be_j = 1;

testrun = zeros(2*Ny,1); % used to compare the order in which the eigenvalues are sorted
for i = 1:2*Ny
    if  ((0.9999 > abs(sort_me(i,i))) || (abs(sort_me(i,i)) > 1.0001))% evanescent
        if (abs(sort_me(i,i)) < 1) % forward
            %s_fe(fe_j,1) = sort_me(i,i);
            %v_fe(:,fe_j) = Vector(:,i); 
            s_fe(1,fe_j) = sort_me(i,i);
            s_fe(2,fe_j) = i;
            fe_j = fe_j+1;
            testrun(i,1) = 3;
        else % backwards
            %s_be(be_j,1) = sort_me(i,i);
            %v_be(:,be_j) = Vector(:,i); 
            s_be(1,be_j) = sort_me(i,i);
            s_be(2,be_j) = i;
            be_j = be_j+1;
            testrun(i,1) = 4;
        end
    else  %propagating
        if (imag(sort_me(i,i)) > 0) % forward
            %s_fp(fp_j,1) = sort_me(i,i);
            %v_fp(:,fp_j) = Vector(:,i);
            s_fp(1,fp_j) = sort_me(i,i);
            s_fp(2,fp_j) = i;
            fp_j = fp_j+1; %%col
            testrun(i,1) = 1;
        else % backward
            %s_bp(bp_j,1) = sort_me(i);
            %v_bp(:,bp_j) = Vector(:,i); 
            s_bp(1,bp_j) = sort_me(i,i);
            s_bp(2,bp_j) = i;
            bp_j = bp_j+1;
            testrun(i,1) = 2;
        end
    end
end

v_fp = Vector(:,double(s_fp(2,:)));
v_fe = Vector(:,double(s_fe(2,:)));
v_bp = Vector(:,double(s_bp(2,:)));
v_be = Vector(:,double(s_be(2,:)));

Vsort = [v_fp v_fe v_bp v_be];
%Dsort = [s_fp(1,:);s_fe(1,:);s_bp(1,:);s_be(1,:)];
Dsort = 0;
fp_modes = s_fp;

end