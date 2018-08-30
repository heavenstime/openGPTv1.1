tsLFileName = '(Directory to MNIST data)/test/t10k-labels-idx1-ubyte';
tsDFileName = '(Directory to MNIST data)/test/t10k-images-idx3-ubyte';
trLFileName = '(Directory to MNIST data)/train/train-labels-idx1-ubyte';
trDFileName = '(Directory to MNIST data)/train/train-images-idx3-ubyte';

outDir      = './gptWork/mnistPgm';

fp = fopen(tsLFileName);
[mgn   cnt] = fread(fp, 1, 'int32', "ieee-be")
[nData cnt] = fread(fp, 1, 'int32', "ieee-be")
[label cnt] = fread(fp, nData, 'uint8');
if cnt != nData
  printf("Data number error\n");
endif

nTs = zeros(1,10);
for cat= 0:9
 nTs(cat + 1) = sum(label == cat);
 printf("%d : %d \n", cat, nTs(cat + 1));
endfor
sum(nTs)

fp = fopen(tsDFileName);
[mgn     cnt] = fread(fp, 1, 'int32', "ieee-be")
[nData   cnt] = fread(fp, 1, 'int32', "ieee-be")
[rows    cnt] = fread(fp, 1, 'int32', "ieee-be")
[columns cnt] = fread(fp, 1, 'int32', "ieee-be")
[data    cnt] = fread(fp, nData * rows * columns, 'uint8');
if cnt != nData * rows * columns
  printf("Data number error\n");
endif

imgData = reshape(data, rows * columns, nData); 
clear data
for cat= 0:9
  catData = imgData(:, find(label == cat));
  for dataN = 1:nTs(cat + 1)
    resData = reshape(255 - catData(:, dataN), rows, columns);
    fileName = sprintf("%s/tdg%1d_%04d_gray.pgm", outDir, cat, dataN);
    writePgm(fileName, resData, rows, columns);
    clear resData
  endfor
  clear catData
endfor  

clear label imgData

% ----------------------------------------------------
fp = fopen(trLFileName);
[mgn   cnt] = fread(fp, 1, 'int32', "ieee-be")
[nData cnt] = fread(fp, 1, 'int32', "ieee-be")
[label cnt] = fread(fp, nData, 'uint8');
if cnt != nData
  printf("Data number error\n");
endif

nTr = zeros(1,10);
for cat= 0:9
 nTr(cat + 1) = sum(label == cat);
 printf("%d : %d \n", cat, nTr(cat + 1));
endfor
sum(nTr)

fp = fopen(trDFileName);
[mgn     cnt] = fread(fp, 1, 'int32', "ieee-be")
[nData   cnt] = fread(fp, 1, 'int32', "ieee-be")
[rows    cnt] = fread(fp, 1, 'int32', "ieee-be")
[columns cnt] = fread(fp, 1, 'int32', "ieee-be")
[data    cnt] = fread(fp, nData * rows * columns, 'uint8');
if cnt != nData * rows * columns
  printf("Data number error\n");
endif

imgData = reshape(data, rows * columns, nData); 
clear data
for cat= 0:9
  catData = imgData(:, find(label == cat));
  for dataN = 1:nTr(cat + 1)
    resData = reshape(255 - catData(:, dataN), rows, columns);
    fileName = sprintf("%s/ldg%1d_%04d_gray.pgm", outDir, cat, dataN);
    writePgm(fileName, resData, rows, columns);
    clear resData
  endfor
  clear catData
endfor  

clear label imgData

