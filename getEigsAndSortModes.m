function [Vsort,Dsort, fp_modes] = getEigsAndSortModes(T0)
    global Ny;
    
    [V, D] = eig(T0);
    
    sort_me = diag(D); %turn diagonal into column
sort_me_real = real(sort_me);
sort_me_imag = imag(sort_me);

s_fp = zeros(2*Ny,1);
s_fe = zeros(2*Ny,1);
s_bp = zeros(2*Ny,1);
s_be = zeros(2*Ny,1);

v_fp = zeros(2*Ny,2*Ny);
v_fe = zeros(2*Ny,2*Ny);
v_bp = zeros(2*Ny,2*Ny);
v_be = zeros(2*Ny,2*Ny);

fp_j = 1;
fe_j = 1;
bp_j = 1;
be_j = 1;


for i = 1:2*Ny
    if ((0.99 < abs(sort_me(i))) && (abs(sort_me(i)) < 1.01)) %propagating
        if (imag(sort_me(i)) > 0) % forward
            s_fp(fp_j,1) = sort_me(i);
            v_fp(:,fp_j) = V(:,i);
            fp_j = fp_j+1;
        else % backward
            s_bp(bp_j,1) = sort_me(i);
            v_bp(:,bp_j) = V(:,i); 
            bp_j = bp_j+1;
        end
    else % evanescent
        if (abs(sort_me(i)) < 1) % forward
            s_fe(fe_j,1) = sort_me(i);
            v_fe(:,fe_j) = V(:,i); 
            fe_j = fe_j+1;
        else % backwards
            s_be(be_j,1) = sort_me(i);
            v_be(:,be_j) = V(:,i); 
            be_j = be_j+1;
        end
    end
end


%for i = 1:2*Ny
%    if ((sort_me_imag(i) > 0)) % x + iy forward propagating
%        s_fp(fp_j,1) = sort_me_real(i) + sort_me_imag(i)*1i;
%        v_fp(:,fp_j) = V(:,i); 
%        fp_j = fp_j+1;
%    elseif (abs(sort_me_real(i) < 1) && (sort_me_imag(i) == 0)) % forward evanescent
%        s_fe(fe_j,1) = sort_me_real(i);
%        v_fe(:,fe_j) = V(:,i); 
%        fe_j = fe_j+1;
%    elseif ((sort_me_imag(i) < 0)) % backward propagating
%        s_bp(bp_j,1) = sort_me_real(i) + sort_me_imag(i)*1i;
%        v_bp(:,bp_j) = V(:,i); 
%        bp_j = bp_j+1;
%    elseif (abs(sort_me_real(i) > 1) && (sort_me_imag(i) == 0)) % backward evanescent 
%        s_be(be_j,1) = sort_me_real(i);
%        v_be(:,be_j) = V(:,i); 
%        be_j = be_j+1;
%    else
%        % Do nothing
%    end
%
%end

zeroRows_fp = any(s_fp==0, 2);
s_fp(zeroRows_fp, :) = [];
v_fp(:,zeroRows_fp) = [];

zeroRows_fe = any(s_fe==0, 2);
s_fe(zeroRows_fe, :) = [];
v_fe(:,zeroRows_fe) = [];

zeroRows_bp = any(s_bp==0, 2);
s_bp(zeroRows_bp, :) = [];
v_bp(:,zeroRows_bp) = [];

zeroRows_be = any(s_be==0, 2);
s_be(zeroRows_be, :) = [];
v_be(:,zeroRows_be) = [];

%zeroRows_Vfp = any(v_fp==0, 2);
%v_fp(zeroRows_Vfp, :) = [];
%zeroRows_Vfe = any(v_fe==0, 2);
%v_fe(zeroRows_Vfe, :) = [];
%zeroRows_Vbp = any(v_bp==0, 2);
%v_bp(zeroRows_Vbp, :) = [];
%zeroRows_Vbe = any(v_be==0, 2);
%v_be(zeroRows_Vbe, :) = [];

%s_fp = sort(s_fp);
%s_fe = sort(s_fe);
%s_bp = sort(s_bp);
%s_be = sort(s_be);
%Vsort = transpose([v_fp;v_fe;v_bp;v_be]);
Vsort = ([v_fp v_fe v_bp v_be]);
Dsort = [s_fp;s_fe;s_bp;s_be];
fp_modes = s_fp;

%out = diag(s);

end