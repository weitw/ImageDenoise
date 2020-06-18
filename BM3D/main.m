clear;clc;
 
pauseTime = 1;

data_path = "..\Set12";
ext = ["*.jpg", "*.png", "*.jpeg"];
filePaths   =  [];
for i = 1 : length(ext)
    filePaths = cat(1,filePaths, dir(fullfile(data_path,ext(i))));
end

noise_leval = [10,15,20,25,30,35,40,45,50,55,60,65,70];

for i = 1:length(noise_leval)
    PSNRs = [];
    SSIMs = [];
    sigma = noise_leval(i);
    for j = 1:length(filePaths)
        y = imread(filePaths(j).name);
        if length(size(y)) > 2
            y = rgb2gray(y);
        end
        y = im2double(y);
        z = y + (sigma/255)*randn(size(y));
        % …˙≥…‘Î…˘ÕºœÒ
        
        [PSNR,SSIM,y_est] = BM3D(y, z, sigma, 'np', 0);
        PSNRs(j) = PSNR;
        SSIMs(j) = SSIM;
    
        imshow(cat(2,im2uint8(y),im2uint8(z),im2uint8(y_est)));
        title([num2str(sigma),'   ', filePaths(j).name,'    ',num2str(PSNR,'%2.2f'),'dB','    ',num2str(SSIMs(j),'%2.4f')])
        drawnow;
        pause(pauseTime)
    end
    disp(["sigma:", sigma, " psnr:", mean(PSNRs), "  ssim:", mean(SSIMs)]);
end



