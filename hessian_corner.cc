/*-------------------------------------------------------------------*/
/*            �w�b�Z�s���Ƃ��̌ŗL�l�v�Z�T���v���v���O����           */
/*                            by H.Nagahashi                         */
/*   �����F���������F���͉摜ARY2I�^�I�u�W�F�N�g                     */
/*         ���������F�ŗL�l�i���j�i�[�p�����^�z���I�u�W�F�N�g        */
/*         ���O�����F�ŗL�l�i���j�i�[�p�����^�z���I�u�W�F�N�g        */
/*         ���l�����F�K�E�X�j�I�u�W�F�N�g                            */
/*   �K�E�X�j�̃X�P�[���p�����[�^��1.0�`3.5���x���I���B              */
/*-------------------------------------------------------------------*/
#include "iip.h"
#include "corner.h"
/*-------------------------------------------------------------------*/
/*                Calculate EigenValue of Hessian                    */
/*-------------------------------------------------------------------*/
void calcHessianEigenVal(ARY2I& img, ARY2D& ary0, ARY2D& ary1, GSKNL0& gk)
{
  int width     = img.width,   height = img.height;
  int2D pic     = img.array;
  double2D eig0 = ary0.array,    eig1   = ary1.array;

  //���Ɨp�����^2�����z���̊m��
  double2D gx2  = allocDouble2D(height,width);
  double2D gx1  = allocDouble2D(height,width);
  double2D gx0  = allocDouble2D(height,width);

  int sX = gk.step[0], sY = gk.step[1]; //�e�����̃}�X�N���Ԉʒu
  ARY2D wx0   = gk.weight[0];  //�������O���A�ꎟ�A�񎟓��֐��p���`�t�B���^�z��
  double2D wx = wx0.array;
  ARY2D wy0   = gk.weight[1];  //�������O���A�ꎟ�A�񎟓��֐��p���`�t�B���^�z��
  double2D wy = wy0.array;

  double sum1,sum2,sum3;
  for(int i=0 ; i < height; i++) {
    for(int j=0; j < width; j++) {
      sum1 = sum2 = sum3 = 0.0;
      for(int k = -sX; k <= sX; k++){
        if((j+k) < 0 || (j+k) >= width) continue;
        sum1 += (double)pic[i][j+k] * wx[0][sX+k];  //���̂O�����֐����ݍ���
        sum2 += (double)pic[i][j+k] * wx[1][sX+k];  //���̈ꎟ���֐����ݍ���
        sum3 += (double)pic[i][j+k] * wx[2][sX+k];  //���̓񎟓��֐����ݍ���
      }
      gx0[i][j] = sum1;
      gx1[i][j] = sum2;
      gx2[i][j] = sum3;
    }
  }
  //�c�����K�E�X�j�ɂ����c�����ݍ��ݏ����ƌŗL�l�v�Z
  double a, b, c;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      a = b = c  = 0.0;
      for(int k = -sY; k <= sY; k++){
	 if((i+k) < 0 || (i+k) >= height) continue;
	 a += gx2[i+k][j] * wy[0][sY+k];
	 b += gx0[i+k][j] * wy[2][sY+k];
         c += gx1[i+k][j] * wy[1][sY+k];
      }
      //�_(i,j)�̃w�b�Z�s�� H=[{a  c}{c  b}]
      //�w�b�Z�s���̌ŗL�l�v�Z ���������� |H-rI|=0
      //r�̂Q�����򎮁Fr^2-(a+b)r+(ab-c^2)=0������
      double dx = (double)(a + b);
      double dy = sqrt((double)((a-b)*(a-b) + c*c) );
      double r0 = (dx + dy )/2.0;           //�傫���ŗL�l
      double r1 = (dx - dy )/2.0;           //�������ŗL�l
      eig0[i][j] = r0;                      //�傫���ŗL�l���i�[
      eig1[i][j] = r1;                      //�������ŗL�l���i�[
    }
  }
  //���Ɨ̈��̉���
  freeDouble2D(gx0, height, width);
  freeDouble2D(gx1, height, width);
  freeDouble2D(gx2, height, width);
}
/*-------------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  //-------------�ۑ�1---------------
  double sigma, gamma;
  if (argc != 5){
    std::cout << "Usage: " << argv[0] << " input output1 output2 sigma gamma" << std::endl;
    return(1);
  }
  PGM& input = readPGM(argv[1]);
  sscanf(argv[3], "%lf", &sigma);
  sscanf(argv[4], "%lf", &gamma);
  printf("%d\n", __LINE__);
  int height = input.body.height, width = input.body.width;
  ARY2D eig0 = createARY2D(height, width);
  ARY2D eig1 = createARY2D(height, width);
  GSKNL0& gk = createGSKNL0(sigma, gamma);
  calcHessianEigenVal(input.body, eig0, eig1, gk);
  //---------------------------------
  //-------------�ۑ�3---------------
  PGM& output = createPGM(height, width);
  extractHessianCorner(output.body, eig0, eig1);
  writePGM(output, argv[2]);
  //---------------------------------
  deletePGM(output);
  deleteARY2D(eig0);
  deleteARY2D(eig1);
  deleteGSKNL0(gk);
  return 0;
}
