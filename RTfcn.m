% RTfcn.m traces rays from a given source (sourcecoords) through a given
% scene (objset) to a given CCD (CCDrange). Outputs captured CCD
% intensities.
function intensity=RTfcn(sourcecoords,CCDrange,objset)
hold off
clc
bounds=[0 −1 2 1]; % scene bounds − [xmin, ymin, xmax, ymax].
sceneon=1; % Plot scene (rays, objects, CCD).
reference=0; % Use reference beam?
maxrefs=12; % Maximum number of reflections + 1.
lamda=5.4e−2; % Wavelength of monochromatic light source.
nrays=800; % Number of rays used.
CCDres=100; % Number of pixels along CCD.

% // Calculate gradients & normal vectors for each surface object
objgrad=(objset(:,4)−objset(:,2))./(objset(:,3)−objset(:,1));
objyint=objset(:,2)−objset(:,1).∗objgrad;
v2=[objset(:,3)−objset(:,1) objset(:,4)−objset(:,2)];
v2n=zeros(size(v2));
for i = 1:size(v2,1)
  v2(i,:) = v2(i,:)/norm(v2(i,:));
  if v2(i,2)==0
    if v2(i,1)==1
      v2n(i,:)=[0 1];
    else
      v2n(i,:)=[0 −1];
    end
  elseif v2(i,1)==0
    if v2(i,2)==−1
      v2n(i,:)=[1 0]; 
    else 
      v2n(i,:)=[−1 0]; 
    end 
  else
    v2n(i,:)=[1,−v2(i,1)/v2(i,2)]; 
    v2n(i,:)=v2n(i,:)/norm(v2n(i,:)); 
  end 
end 

%// Setup ray matrix: 
% Column key: (1)normalvect[x],(2)normalvect[y],(3)gradient,(4)startcoords[x],(5)startcoords[y], 
% (6)endcoords[x],(7)endcoords[y],(8)obj.int,(9)yint 
raymat=zeros(nrays,9,maxrefs); 
if sourcecoords(4)==sourcecoords(2) % Special case when beam points vertically down. 
  raymat(:,1,1)=0; 
  raymat(:,2,1)=−1; 
  raymat(:,3,1)=−inf; 
  raymat(:,4,1)=linspace(sourcecoords(1),sourcecoords(3),nrays); 
  raymat(:,5,1)=sourcecoords(2); 
  raymat(:,6,1)=raymat(:,4,1); 
  raymat(:,7,1)=bounds(2); 
else
  g1=−1/((sourcecoords(4)−sourcecoords(2))/(sourcecoords(3)−sourcecoords(1))); 
  raymat(:,1:2,1)=repmat([1 g1]/norm([1 g1]),nrays,1); 
  raymat(:,3,1)=g1;
  raymat(:,4,1)=linspace(sourcecoords(1),sourcecoords(3),nrays);
  raymat(:,5,1)=linspace(sourcecoords(2),sourcecoords(4),nrays);
  raymat(:,9,1)=raymat(:,5,1)−raymat(:,3,1).∗raymat(:,4,1);
  kkk = raymat;
end 

%// Update raymatrix according to layout of the scene: 
raymat=computerefs(raymat,objset,objgrad,objyint,v2n,bounds); 

%// REFERENCE BEAM 
if reference==1 
  ref.coords=[sourcecoords(3),sourcecoords(2),sourcecoords(1),sourcecoords(4)]; 
  ref.grad = −1/((ref.coords(4)−ref.coords(2))/(ref.coords(3)−ref.coords(1))); 
  ref.nrays=nrays; 
  refmat=zeros(ref.nrays,9,maxrefs); 
  refmat(:,1:2,1)=repmat([1 ref.grad]/norm([1 ref.grad]),ref.nrays,1); 
  refmat(:,3,1)=ref.grad; 
  refmat(:,4,1)=linspace(ref.coords(1),ref.coords(3),ref.nrays); 
  refmat(:,5,1)=linspace(ref.coords(2),ref.coords(4),ref.nrays); 
  refmat(:,9,1)=refmat(:,5,1)−refmat(:,3,1).∗refmat(:,4,1); 
  refmat(:,6,1)=bounds(3); 
  refmat(:,7,1)=refmat(:,3,1)∗bounds(3)+refmat(:,9,1); 
  CCDrefrays=cell(CCDres,4); 
  refmat=computerefs(refmat,objset,objgrad,objyint,v2n,bounds); 
end 

%// Draw objects, sources, rays & CCD (can be very time consuming). 
if sceneon==1 
  for i = 1:size(objset,1); 
    line([objset(i,1),objset(i,3)],[objset(i,2),objset(i,4)],"color","k");
    hold on;
  end
for i=1:nrays
  for k=1:maxrefs−1 
    line([raymat(i,4,k),raymat(i,6,k)],[raymat(i,5,k),raymat(i,7,k)]) 
  end 
end 
line([sourcecoords(1),sourcecoords(3)],[sourcecoords(2),sourcecoords(4)],"color","k"); 
line([CCDrange(1) CCDrange(3)],[CCDrange(2),CCDrange(4)],"color","k") 
if reference==1 
  line([ref.coords(1),ref.coords(3)],[ref.coords(2),ref.coords(4)],"color","k"); 
    for i=1:ref.nrays 
      for k = 1:maxrefs−1 
        line([refmat(i,4,k),refmat(i,6,k)],[refmat(i,5,k),refmat(i,7,k)],"color","r"); 
      end 
    end 
  end 
axis([bounds(1) bounds(3) bounds(2) bounds(4)]) 
end 

%// CAMERA 
CCDrays=cell(CCDres,4); % For each pixes, contains: [which rays hit it, which parts of each ray hit it, OPL 
% of each ray hitting it, phase of each ray hitting it]. 
inc=(CCDrange(4)−CCDrange(2))/CCDres; 
if CCDrange(3)˜=CCDrange(1) 
  inc = sqrt((CCDrange(3)−CCDrange(1))ˆ2+(CCDrange(4)−CCDrange(2))ˆ2)/CCDres; 
  CCDgrad=(CCDrange(4)−CCDrange(2))/(CCDrange(3)−CCDrange(1)); 
  CCDyint=CCDrange(2)−CCDgrad∗CCDrange(1);
end
for j=1:CCDres
  xmin=min(CCDrange(1)+(j−1)∗(CCDrange(3)−CCDrange(1))/CCDres,CCDrange(1)+j∗(CCDrange(3)−CCDrange(1))/CCDres);
  xmax=max(CCDrange(1)+(j−1)∗(CCDrange(3)−CCDrange(1))/CCDres,CCDrange(1)+j∗(CCDrange(3)−CCDrange(1))/CCDres); 
  ymin=min(CCDrange(2)+(j−1)∗(CCDrange(4)−CCDrange(2))/CCDres,CCDrange(4)+j∗(CCDrange(4)−CCDrange(2))/CCDres); 
  ymax=max(CCDrange(2)+(j−1)∗(CCDrange(4)−CCDrange(2))/CCDres,CCDrange(4)+j∗(CCDrange(4)−CCDrange(2))/CCDres); 
  for k = 1:size(raymat,3)−1 % Compute hits on each pixel. 
    for i = 1:size(raymat,1) 
      v=raymat(i,3,k)∗CCDrange(1)+raymat(i,9,k); 
      if CCDrange(3)==CCDrange(1) && v>=CCDrange(2)+(j−1)∗inc && v<= CCDrange(2)+j∗inc && min(raymat(i,4,k),raymat(i,6,k))<=CCDrange(1) && max(raymat(i,4,k),raymat(i,6,k))>=CCDrange(3) 
        if k==1 | (k>1 & raymat(i,6:7,k)˜=raymat(i,6:7,k−1)) 
          CCDrays{j,1}=[CCDrays{j,1},i]; 
          CCDrays{j,2}=[CCDrays{j,2},k]; 
          CCDrays{j,3}=[CCDrays{j,3},OPL(raymat(i,:,:),k,CCDrange(1))]; 
          CCDrays{j,4}=[CCDrays{j,4},mod(OPL(raymat(i,:,:),k,CCDrange(1)),lamda)∗2∗pi/lamda]; 
        end 
      elseif CCDrange(3)˜= CCDrange(1) 
        CCDint=(CCDyint−raymat(i,9,k))/(raymat(i,3,k)−CCDgrad); 
        % Below: Does the ray hit the pixel? EDIT add a CCDyint to 
        % simplify. 
        if CCDint>=xmin && CCDint<=xmax && raymat(i,3,k)∗CCDint+raymat(i,9,k)>=ymin && raymat(i,3,k)∗CCDint+raymat(i,9,k)<=ymax && CCDint>=min(raymat(i,4,k),raymat(i,6,k)) && CCDint <= max(raymat(i,4,k),raymat(i,6,k)) && CCDint∗raymat(i,3,k)+raymat(i,9,k)>=min(raymat(i,5,k),raymat(i,7,k)) && CCDint∗raymat(i,3,k)+raymat(i,9,k)<=max(raymat(i,5,k),raymat(i,7,k)) 
          if k==1 | (k>1 & raymat(i,6:7,k)˜=raymat(i,6:7,k−1)) % Prevents rays stopped at the boundary from solving 
            CCDrays{j,1}=[CCDrays{j,1},i]; % repeatedly with remaining reflection iterates.
            CCDrays{j,2}=[CCDrays{j,2},k];
            CCDrays{j,3}=[CCDrays{j,3},OPL(raymat(i,:,:),k,CCDint)];
            CCDrays{j,4}=[CCDrays{j,4},mod(OPL(raymat(i,:,:),k,CCDint),lamda)∗2∗pi/lamda];
        end
      end 
    end 
  end 
end 
if reference==1 
  for k = 1:(size(refmat,3)−1) 
    for i = 1:size(refmat,1) 
      v=refmat(i,3,k)∗CCDrange(1)+refmat(i,9,k); 
        if CCDrange(3)==CCDrange(1) && v>=CCDrange(2)+(j−1)∗inc && v< CCDrange(2)+j∗inc && min(refmat(i,4,k), refmat(i,6,k))<=CCDrange(1) && max(refmat(i,4,k),refmat(i,6,k))>=CCDrange(3) 
          if k==1 | (k>1 & refmat(i,6:7,k)˜=refmat(i,6:7,k−1)) 
            CCDrefrays{j,1}=[CCDrefrays{j,1},i]; 
            CCDrefrays{j,2}=[CCDrefrays{j,2},k]; 
            CCDrefrays{j,3}=[CCDrefrays{j,3},OPL(refmat(i,:,:),k,CCDrange(1))]; 
            CCDrefrays{j,4}=[CCDrefrays{j,4},mod(OPL(refmat(i,:,:),k,CCDrange(1)),lamda)∗2∗pi/lamda]; 
          end 
        elseif CCDrange(3)˜= CCDrange(1) 
          CCDint=(CCDyint−refmat(i,9,k))/(refmat(i,3,k)−CCDgrad); 
          if CCDint>=xmin && CCDint<=xmax && refmat(i,3,k)∗CCDint+refmat(i,9,k)>=ymin && refmat(i,3,k)∗CCDin+refmat(i,9,k)<=ymax && CCDint>=min(refmat(i,4,k),refmat(i,6,k)) && CCDint <= max(refmat(i,4,k),refmat(i,6,k)) && CCDint∗refmat(i,3,k)+refmat(i,9,k)>=min(refmat(i,5,k),refmat(i,7,k)) && CCDint∗refmat(i,3,k)+refmat(i,9,k)<=max(refmat(i,5,k),refmat(i,7,k))
            if k==1 | k>1 && refmat(i,6:7,k)˜=refmat(i,6:7,k−1)
              CCDrefrays{j,1}=[CCDrefrays{j,1},i];
              CCDrefrays{j,2}=[CCDrefrays{j,2},k];
              CCDrefrays{j,3}=[CCDrefrays{j,3},OPL(refmat(i,:,:),k,CCDint)]                  CCDrefrays{j,4}=[CCDrefrays{j,4},mod(OPL(refmat(i,:,:),k,CCDint),lamda)∗2∗pi/lamda]; 
            end 
          end 
        end 
      end 
    end 
  end    
end 
intensity=zeros(CCDres,1); % Compute intensity of each pixel. 
if reference==1     
  for i = 1:CCDres 
    intensity(i,1)=abs(sum(exp(1i∗opl2phase(CCDrays{i,3},lamda)))+sum(exp(1i∗opl2phase(CCDrefrays{i,4},lamda))));
  end 
else
  for i = 1:CCDres          
    intensity(i,1)=abs(sum(exp(1i∗opl2phase(CCDrays{i,3},lamda)))); 
  end 
end 
if CCDon==1 
  figure 
  imshow(repmat(flipud(intensity),1,CCDres))     
  caxis([0 5])
end
