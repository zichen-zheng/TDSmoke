function gray2bit(imgpath, outpath)

img = imread(imgpath);
out = img > 127;
fout = fopen(outpath, 'w');
for i = 1 : size(img, 1)
    for j = 1 : size(img, 2)
        fprintf(fout, '%d ', out(i,j));
    end;
    fprintf(fout, '\n');
end
