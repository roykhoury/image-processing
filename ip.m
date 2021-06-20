% Variables
weak = 100;
strong = 255;

% Read and transform image to grayscale
rgb = imread('test.png');
gray = rgb2gray(rgb);

% Apply a gaussion filter on the image in order to reduce noise
gaussian = gaussBlur(gray, size(gray, 1), size(gray, 2), 2);

% Apply a gradient Sobel Filter on the image to start detecting edges
[gradient, direction] = sobelFilter(gaussian);

% Thin out the edges using Non-Maximum Suppression
suppressed = nonMaxSuppression(gradient, direction);

% Change all lines in the image to either strong or weak lines
threshold = doubleThreshold(suppressed, 0.1, 0.2, weak, strong);

% Identify relevant weak edges and add them to image
final = weakEdgeTracking(threshold, weak, strong);

imshow(final)

%% Gaussian Filter
% Reduce image noise by applying a Gaussian filter 
% - source channel (the image source matrix)
% - filter width
% - filter height 
% - StDev (radius or sigma)
% This algorithm is slower then the others, however it perfectly applies
% the filter with average error per pixel of 0.00%
function out = gaussBlur(scl, w, h, r)
    rs = ceil(r * 2.57);
    for i = 1:h
        for j = 1:w
            val = 0;
            wsum = 0;
            for iy = i-rs:i+rs+1
                for ix = j-rs:j+rs+1
                    x = min(w-1, max(0, ix));
                    y = min(h-1, max(0, iy));
                    dsq = (ix-j)*(ix-j)+(iy-i)*(iy-i);
                    wght = exp( -dsq / (2*r*r) ) / (pi*2*r*r);
                    val = val + scl(x+1, y+1) * wght;  
                    wsum = wsum + wght;
                end
            out(j, i) = round(val/wsum);            
            end
        end
    end
end


%% Sobel Filter
% Apply a sobel filter on the pre-blurred image to detect edges and highlight
% them in a white color, while reducing the other pixels to black
% - Image matrix to detect edges for
function [out, angle] = sobelFilter(img)
    Kx = [-1 0 1; -2 0 2; -1 0 1]; Ky = Kx';
    
    Ix = conv2(double(img), Kx, 'same');
    Iy = conv2(double(img), Ky, 'same');
    
    out = uint8(sqrt(Ix.^2 + Iy.^2));
    angle = uint8(atan2(Iy, Ix)*180/pi);
end

%% Non-Maximum Suppression
% Thing out the edges of the edge image
% This is ideally done since the final image should have thin uniform lines
% - The img matrix
% - The edge directions
function out = nonMaxSuppression(img, angle)
    [M, N] = size(img);
    Z = zeros(M,N);
    
    if (angle < 0)
        angle = angle + 180;
    end
    
    for i = 1:(M-1)
        for j = 1:(N-1)
            try
                q = 255;
                r = 255;
                
                % angle 0
                if (0 <= angle(i,j) && angle(i,j) < 22.5) || (157.5 <= angle(i,j) && angle(i, j) <= 180)
                    q = img(i, j+1);
                    r = img(i, j-1);
                % angle 45
                elseif (22.5 <= angle(i,j) && angle(i, j) < 67.5)
                    q = img(i+1, j-1);
                    r = img(i-1, j+1);
                % angle 90
                elseif (67.5 <= angle(i,j) && angle(i, j) < 112.5)
                    q = img(i+1, j);
                    r = img(i-1, j);
                % angle 135
                elseif (112.5 <= angle(i,j) && angle(i, j) < 157.5)
                    q = img(i-1, j-1);
                    r = img(i+1, j+1);
                end
                
                if (img(i,j) >= q) && (img(i,j) >= r)
                    Z(i,j) = img(i,j);
                else
                    Z(i,j) = 0;
                end

            catch E 
            end
        end
    end
    out = uint8(Z);
end

%% Double Threshhold
% Identify the img's strong, weak and non-relevant lines
% - The img matrix
% - The lower end threshold ratio
% - The higher end threshold ratio
function out = doubleThreshold(img, low, high, weak, strong)
    highThreshold = max(max(img)) * high;
    lowThreshold = highThreshold * low;
    
    for i=2:size(img, 1) - 1 % row
        for j=2:size(img, 2) - 1 % col
            if img(i,j) > highThreshold    % Strong edge
                img(i,j) = strong;
            elseif img(i,j) < lowThreshold % No edge
                img(i,j) = 0;
            else                            % Weak edge
                img(i, j) = weak;
            end
        end
    end
    out = uint8(img);
end

%% Weak edge tracking and filtering
% Identify the weak edges that are actually relevant in the image and
% change them into strong edges, otherwise remove them
% - The image source matrix
function out = weakEdgeTracking(img, weak, strong)
    [M, N] = size(img);
    for i = 1:M-1
        for j = 1:N-1
            if (img(i,j) == weak)
                try
                    if ((img(i+1, j-1) == strong) || (img(i+1, j) == strong) || (img(i+1, j+1) == strong) || (img(i, j-1) == strong) || (img(i, j+1) == strong) || (img(i-1, j-1) == strong) || (img(i-1, j) == strong) || (img(i-1, j+1) == strong))
                        img(i, j) = strong;
                    else
                        img(i, j) = 0;
                    end
                catch EX
                end
            end
        end
    end
    out = uint8(img);
end
