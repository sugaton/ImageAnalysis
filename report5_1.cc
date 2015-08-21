#include "iip.h"
#include <cmath>

int main(int argc,char *argv[])
{
int type;
if (argc!=4 ||(argc>1 && strcmp(argv[1],"help")==0)) {
cout << "\nUsage: " << argv[0] << " inImage outImage type\n" <<endl;
exit(1);
}
sscanf(argv[3], "%d", &type);
PGM& pgm = readPGM(argv[1]);
PGM& outpgm = createPGM(pgm.body.height, pgm.body.width);
switch (type) {
  case 1:
    rankMedian(pgm.body, outpgm.body);
  case 2:
    rankMin(pgm.body, outpgm.body);
  case 3:
    rankMax(pgm.body, outpgm.body);
  case 4:
    localAverageARY2I(pgm.body, outpgm.body);
}
writePGM(outpgm, argv[2]);
deletePGM(outpgm);
deletePGM(pgm);
return 0;
}
