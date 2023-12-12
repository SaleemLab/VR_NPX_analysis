classdef LFPanalysis
    %LFPANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        function [CSDoutput] = CSD(LFP,fs,spacing)
            fs=fs/1000;           % sampling frequency (kHz)
            chans=1:size(LFP,2);     % channel used
            el_spac=spacing/1000;    % electrode spacing (mm);
            cond_top=0;     % Conductivity at the surface of the brain (oil ~ 0 S/m)
            
            el_pos=((chans)-1)*el_spac+0.01; % electrode positions in mm from surface
            
            % Set a bunch of parameters
            % Electrode contacts positions along the z axis (mm from surface)
            Par.el_pos=el_pos;
            % Conductivity of the cortical tissue along the z axis: sigmaz (S/m)
            Par.cond=0.3;
            % Conductivity at the surface of the brain (oil ~ 0)  (S/m)
            Par.cond_top=cond_top;
            % Diameter of the cortical column activated by the stimulus (mm)
            Par.diam=0.5;
            % Three Point Filter parameters (Hamming filter)
            Par.b0=0.54; % center
            Par.b1=0.23; % neighbour
            
            csd = LFPanalysis.computeCSD(LFP,Par);
            CSDoutput = csd.data;
        end
        % CSD functions
        function out=computeCSD(LFPs, Par)
            % out=computeCSD(LFPs, Par)
            %
            % This is the Delta-Source iCSD method with spatial smoothing
            %
            % It calls routines extracted from the GUI of the iCSD toolbox
            % available here : http://arken.umb.no/~klaspe/iCSD.php
            %
            % The method and its advantages are described in details in:
            %
            % K. Pettersen, A. Devor, I. Ulbert, A.M. Dale and G.T. Einevoll,
            % Current-source density estimation based on inversion of electrostatic
            % forward solution: Effects of finite extent of neuronal activity and
            % conductivity discontinuities, Journal of Neuroscience Methods (2006).
            %
            % out contains the resulting CSD in uA/mm3
            %
            %  Requirements for the LFPs matrix:
            %   - accepts ANY number of dimensions, as long as...
            %   - the last dimension is depth/channel
            %   - For each site, the mean of the LFP across time is zero.
            %   - LFP is in mV.
            %
            %  Par is a structure containing the necessary parameters:
            %   -b0         Three Point Filter center parameter
            %   -b1         Three Point Filter neighbour parameter
            %   -cond:      cortical conductivity
            %   -cond_top:  conductivity on top of cortex
            %   -diam:      activity diameter
            %   -el_pos:    the z-positions of the electrode contacts
            %
            %Copyright 2010 François D. Szymanski under the General Public License,
            %
            %This program is free software; you can redistribute it and/or
            %modify it under the terms of the GNU General Public License
            %as published by the Free Software Foundation; either version 2
            %of the License, or any later version.
            %
            %See: http://www.gnu.org/copyleft/gpl.html
            
            
            %~~~~~~~~~~~~~~LOAD PARAMS ~~~~~~~~~~~~~~~~~~~~~~
            % Electrode contacts positions along the z axis (mm from surface)
            el_pos=Par.el_pos;
            % Conductivity of the cortical tissue along the z axis (S/m)
            cond=Par.cond;
            % Conductivity at the surface of the brain (S/m)
            cond_top=Par.cond_top;
            % Diameter of the cortical column activated by the stimulus (mm)
            diam=Par.diam;
            % Three Point Filter parameters
            b0=Par.b0; % center
            b1=Par.b1; % neighbour
            %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            
            % Note: LFPs last dimension is depth.
            % However, Pettersen's function accept LFP matrices that are DepthXTime.
            
            % LFPs: 4D with last dim Depth
            % LFPtmp: same data in a 2D matrix with 1st Dim Depth
            LFPtmp=reshape(LFPs,[numel(LFPs)/size(LFPs,ndims(LFPs)) size(LFPs,ndims(LFPs))]);
            LFPtmp=shiftdim(LFPtmp,1);
            
            % We convert from mV to V
            LFPtmp=LFPtmp/1000;
            
            % Compute delta iCSD with Pettersen's routines
            % Results in A/m3.
            CSD=LFPanalysis.deltaCSD(b0,b1,cond,cond_top,diam,el_pos(1:size(LFPtmp,1)),LFPtmp);
            
            % CSD: 2D matrix with 1st Dim Depth
            % transforms to 4D with last dim Depth
            newsize=size(LFPs);
            newsize(end)=size(CSD,1);
            CSD=shiftdim(CSD,1);
            CSD=reshape(CSD,newsize);
            
            %Convert from A/m3 to uA/mm3
            CSD=CSD/1000;
            
            out.data=CSD;
            out.csdpar=Par;
            
        end
        function CSD=deltaCSD(b0,b1,cond,cond_top,diam,el_pos,LFPs)
            %function [CSD,el_pos]=deltaCSD(b0,b1,cond,cond_top,diam,el_pos,LFPs)
            %
            %Computes the filtered delta iCSD
            %b0         Three Point Filter center parameter
            %b1         Three Point Filter neighbour parameter
            %cond:      cortical conductivity, default: 0.3
            %cond_top:  conductivity on top of cortex, default: cond
            %diam:      activity diameter, default: 500e-6
            %el_pos:    the z-positions of the electrode contacts
            %
            %Copyright 2005 Klas H. Pettersen under the General Public License,
            %
            %This program is free software; you can redistribute it and/or
            %modify it under the terms of the GNU General Public License
            %as published by the Free Software Foundation; either version 2
            %of the License, or any later version.
            %
            %See: http://www.gnu.org/copyleft/gpl.html
            
            % filter parameters:
            if b0+2*b1 == 0 & b1~=0
                errordlg('Singularity: b0+2*b1 cannot equal zero.');
                return
            end;
            
            % electrical parameters:
            if cond<=0
                errordlg('Ex. cond. has to be a positive number');
                return
            end;
            
            % size, potential (m1 has to equal number of electrode contacts)
            [m1,m2] = size(LFPs);
            
            % geometrical parameters:
            diam = diam*1e-3; %diameter in [m]
            if diam<=0
                errordlg('Diameter has to be a positive number.');
                return
            end;
            
            el_pos = el_pos*1e-3;
            if cond_top~=cond & (el_pos~=abs(el_pos) | length(el_pos)~=length(nonzeros(el_pos)))
                errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
                return;
            end;
            if m1~=length(el_pos)
                errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
                    num2str(length(el_pos)),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.'])
                return
            end;
            
            % compute delta iCSD:
            CSD = LFPanalysis.F_delta(el_pos,diam,cond,cond_top)^-1*LFPs;
            
            if b1~=0 %filter iCSD
                [n1,n2]=size(CSD);
                CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
                CSD_add(n1+2,:)=zeros(1,n2);
                CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
                CSD = LFPanalysis.S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
            end;
        end
        function out = F_delta(el_pos,d,cond,cond_top)
            %function out = F_delta(el_pos,d,cond,cond_top)
            %
            %Computes the F-matrix from infinitesimally thin current source density
            %sheets with diameter d and homogenous activity throughout the sheet.
            %
            %el_pos:    the z-positions of the electrode contacts, default:
            %100e-6:100e-6:2300e-6
            %d:         activity diameter, default: 500e-6
            %cond:      cortical conductivity, default: 0.3
            %cond_top:  conductivity on top of cortex, default: cond
            %
            %out is a (number_of_electrodes)x(number_of_electrodes) matrix.
            
            %Copyright 2005 Klas H. Pettersen under the General Public License,
            %
            %This program is free software; you can redistribute it and/or
            %modify it under the terms of the GNU General Public License
            %as published by the Free Software Foundation; either version 2
            %of the License, or any later version.
            %
            %See: http://www.gnu.org/copyleft/gpl.html
            
            if nargin < 1, el_pos = 100e-6:100e-6:2300e-6; end
            if nargin < 2, d = 500e-6; end;
            if nargin < 3, cond = 0.3; end;
            if nargin < 4, cond_top = cond; end;
            
            N = length(el_pos);
            z1 = el_pos(1);
            h = el_pos(2)-z1;
            
            try
                load(full_filename,'Fd');
                out = Fd;
            catch
                msgstr = lasterr;
                out = zeros(N);
                for j=1:N                     %zj is position of CSD-plane
                    zj = z1 + (j-1)*h;
                    for i=1:N                   %zi is position of electrode
                        zi = z1 + (i-1)*h;
                        out(j,i) = h/(2*cond)*((sqrt((zj-zi)^2+(d/2)^2)-abs(zj-zi))+ ...
                            (cond-cond_top)/(cond+cond_top)*(sqrt((zj+zi)^2+(d/2)^2)-abs(zj+zi)));
                    end;
                end;
            end;
        end
        function out = S_general(N,b0,b1)
            %S = S_general(N,b0,b1)
            %This is the three point filter matrix.
            %Returns matrix of size (N-2)x(N),
            %which represents a three point "spatial noise" filter with mid
            %coeffescient b0 (will be normalized by the function) and neighbouring
            %coeffescients b1 (will also be normalized).
            %Default filter has b0 = 2 and b1 = 1 and number_of_electrodes = 20.
            %
            %The Hamming-filter has b0 = 0.54 and b1 = 0.23.
            
            %Copyright 2005 Klas H. Pettersen under the General Public License,
            %
            %This program is free software; you can redistribute it and/or
            %modify it under the terms of the GNU General Public License
            %as published by the Free Software Foundation; either version 2
            %of the License, or any later version.
            %
            %See: http://www.gnu.org/copyleft/gpl.html
            
            if nargin < 1, N = 20; end
            if nargin < 3, b1 = 1; b0 = 2; end;
            
            c = b0 + 2*b1;
            
            out = zeros(N-2,N);
            for i=1:N-2
                for j=1:N
                    if (i == j-1)
                        out(i,j) = b0/c;
                    elseif (abs(i-j+1) == 1)
                        out(i,j) = b1/c;
                    end;
                end;
            end;
        end
    end
end