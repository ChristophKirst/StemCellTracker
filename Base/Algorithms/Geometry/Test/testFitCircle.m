%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test fitting circles %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fitting a cricle

% data:
n = 15;
phi = sort(2 * pi * rand(n,1));
r   = 3 + 0.3 * rand(n,1);

x = r .* sin(phi);
y = r .* cos(phi);

figure(1); clf
plot(x,y, '*')


%%

[c,r,e] =fitCircle([x'; y'])

%%

B = [x.^2 + y.^2, x, y, ones(n,1)];
[U,S,V] = svd(B);

se1 = V(:,end)

p = se1(1)
q = se1(2:3)
s = se1(4)


%% plot the result
figure(1); clf
plot(x,y, '*')
viscircles(c',r)





