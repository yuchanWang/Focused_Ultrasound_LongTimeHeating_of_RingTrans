clear;
load('PressureField.mat');
p = press;
target = 0;
for x = 1:490
    for y = 1:490
        if((x-245)^2+(y-245)^2 >= 223^2)
            p(x, y) = 0;
        end
        if p(x,y)~= 0
            p(x, y) = log10(p(x, y));
        end
%         if p(x,y)<target
%             target = p(x, y);
%         end
    end
end
figure;imagesc([0:22:220],[0:22:220],p);axis xy;
xlabel('x position (mm)');ylabel('y position (mm)');


