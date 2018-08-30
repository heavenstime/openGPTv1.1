function writePgm(fileName, val, nx, ny)
  fp = fopen(fileName, 'w');
  fprintf(fp, "P5\n");
  fprintf(fp, "#MNIST %s\n", fileName);
  fprintf(fp, "%d %d\n", nx, ny);
  fprintf(fp, "255\n", fileName);
  
  val1 = (val >= 0) .* val;
  flg2 = (val1 > 255);
  val2 = flg2 * 255 + (1 - flg2) .* val1;
  fwrite(fp, val2, 'uint8');

  fclose(fp);
end