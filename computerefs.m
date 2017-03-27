% computerefs.m updates the raymatrix according to the layout of the scene.
% Additional data such as objgrad and v2n are passed through to save
% computation. bounds is required to prevent rays going to infinity when
% leaving the scene.
function raymat=computerefs(raymat,objset,objgrad,objyint,v2n,bounds) 
for k = 1:size(raymat,3)−1
  hitmat=zeros(size(raymat,1),size(objset,1)); % Will contain x−coordinate of every possible intersection. %// Determine which rays hit which objects
  for n=1:size(raymat,1)
    for m=1:size(objset,1)
      if objset(m,3)==objset(m,1)    % Vertical surfaces (infinite gradient) require special treatment.
        xint = objset(m,1);
        if raymat(n,3,k)∗xint+raymat(n,9,k)>=min(objset(m,2),objset(m,4)) && raymat(n,3,k)∗xint+raymat(n,9,k)<max(objset(m,4),objset(m,2))
          hitmat(n,m)=xint; % if statement requires intersection to occur within the line bounds.
        end
      else % Non−vertical case.
        xint = (−raymat(n,9,k)+(objset(m,2)−objgrad(m)∗objset(m,1)))./(raymat(n,3,k)−objgrad(m));
        yint = raymat(n,3,k)∗xint+raymat(n,9,k);
        if(min(objset(m,1),objset(m,3))<=xint && xint<=max(objset(m,1),objset(m,3))) && yint>=min(objset(m,4),objset(m,2)) && yint <=max(objset(m,4),objset(m,2)); 
          hitmat(n,m)=xint;
        end 
      end
    if abs(raymat(n,3,k))==inf && raymat(n,4,k)>=min(objset(m,1),objset(m,3)) && raymat(n,4,k)<=max(objset(m,1),objset(m,3))
      hitmat(n,m)=raymat(n,4,k); % Vertical rays require special treatment.
    end
    if sign(hitmat(n,m)−raymat(n,4,k))~=sign(raymat(n,1,k)) && hitmat(n,m)−raymat(n,4,k)~=0 
      hitmat(n,m)=0; % Clear entry if ray is travelling away from the surface.
%// Find which object, if any, is hit first for each ray & calculate appropriate reflectance gradients.
HM=zeros(size(hitmat));
for i=1:size(hitmat,1)
  for j=1:size(hitmat,2)
    if hitmat(i,j)~=0
  end 
end
HM(abs(HM)<1e−12)=0; 
HM(HM<=0)=inf;
% Removes the starting point from possible hits, incase of numerical error.
% This will allow the closest hit to be found.

for i=1:size(raymat,1)
  [Y I]=min(HM(i,:));
  if Y~=inf
    raymat(i,8,k)=I;
    if k>1 && raymat(i,8,k)==raymat(i,8,k−1)
     raymat(i,8,k)=0;    % Needed to stop rays at bounds.
    end
  end
  if raymat(i,8,k)~=0
    if raymat(i,3,k)~=−inf
      raymat(i,6,k)=hitmat(i,I); % Terminate ray at hitpoint. 
      raymat(i,7,k)=raymat(i,3,k)∗raymat(i,6,k)+raymat(i,9,k);
    else 
      raymat(i,6,k)=raymat(i,4,k); 
      raymat(i,7,k)=objgrad(raymat(i,8,k))∗raymat(i,4,k)+objyint(raymat(i,8,k));
    end
  if abs(raymat(i,6,k))==inf     % Terminate ray at bounds if need be.
    if sign(raymat(i,1,k))==1
      raymat(i,7,k)=bounds(4);
    else
      raymat(i,7,k)=bounds(2);
    end 
  end
else
  raymat(i,6,k)=bounds(3);
  raymat(i,7,k)=bounds(3)∗raymat(i,3,k)+raymat(i,9,k); else
  raymat(i,6,k)=bounds(1);
  raymat(i,7,k)=bounds(1)∗raymat(i,3,k)+raymat(i,9,k); end
end
if abs(raymat(i,7,k))==inf
  if sign(raymat(i,2,k))==1
    raymat(i,7,k)=bounds(4);
    raymat(i,7,k)=bounds(2);
%// Compute normal vector, gradient, start points and % y−intercept for reflected ray.
  end 
end

if k>1 && (raymat(i,6,k−1)==bounds ( 1 ) | | raymat(i,6,k−1)==bounds(3) | | raymat(i,5,k−1)==bounds(4))
  raymat(i,1:2,k+1)=−2∗dot(raymat(i,[1 2],k),v2n(raymat(i,8,k),:))∗v2n(raymat(i,8,k),:)+raymat(i,[1 2],k); 
  raymat(i,3,k+1)=raymat(i,2,k+1)/raymat(i,1,k+1);
  raymat(i,4:5,k+1)=raymat(i,6:7,k);
  raymat(i,9,k+1)=raymat(i,5,k+1)−raymat(i,3,k+1)∗raymat(i,4,k+1);
if sign(raymat(i,1,k+1))==1 % Ensure ray begins by terminating at bounds (hits computed on next run). raymat(i,6,k+1)=bounds(3);
  raymat(i,7,k+1)=raymat(i,3,k+1)∗bounds(3)+raymat(i,9,k+1);
elseif sign(raymat(i,1,k+1))==−1
  raymat(i,6,k+1)=bounds(1); 
  raymat(i,7,k+1)=raymat(i,3,k+1)∗bounds(1)+raymat(i,9,k+1);
else
% EDIT may need y−direction conditions. 
  raymat(i,6,k+1)=raymat(i,4,k+1); raymat(i,7,k+1)=bounds(4);
  raymat(i,5,k−1)==bounds(2) 
  raymat(i,6:7,k)=raymat(i,4:5,k); % Ensures when a ray terminates at bounds it stays that way.
  end
  if k == 1 && raymat(i,8,k)==0
    raymat(i,7,k)=bounds(2); 
    if raymat(i,3,k)==−inf
      raymat(i,6,k)=raymat(i,4,k);
    else
      raymat(i,6,k)=(raymat(i,7,k)−raymat(i,9,k))/raymat(i,3,k); 
    end
  end
end
