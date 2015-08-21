#include "iip.h"
/*-------------------------------------------------------------------*/
/*   �w�b�Z�s��̌ŗL�l����R�[�i�[�_�����o����T���v���v���O����    */
/*                       by H.Nagahashi                              */
/*      �������FPGM�摜�I�u�W�F�N�g                                */
/*      �������F�ŗL�l�i��j�z��I�u�W�F�N�g                       */
/*      ��O�����F�ŗL�l�i���j�z��I�u�W�F�N�g                       */
/*-------------------------------------------------------------------*/
inline bool is_edge(double val0, double val1, double thre1, double thre2)
{
    return (val0 < thre1) && (val0 >= thre2) && (val1 >= thre2);
}

void extractHessianCorner(ARY2I& img, ARY2D& ary0, ARY2D& ary1)
{
  int height    = img.height,  width = img.width;
  int2D pic     = img.array;
  double2D eig0 = ary0.array,   eig1 = ary1.array;

  //�Ǐ��œK����T�����邽�߂̂Q�����z��
  int2D tmp = allocInt2D(height, width);

  //�傫���ŗL�l�̍ő�A�ŏ������o
  double mn=10000.0, mx = -10000.0;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      if(eig0[i][j] > mx) mx = eig0[i][j];
      if(eig0[i][j] < mn) mn = eig0[i][j];
    }
  }
  // �傫���ŗL�l�̍ő�l��60��(thresh1)�ȏ�����_���G�b�W�����Ƃ��ċp��
  // �傫���ŗL�l�̍ő�l��5��(thresh2)������臒l�Ƃ��Č���
  // �傫���ŗL�l�Ə������ŗL�l��������thresh2�ȏ�ł����f���p�_�Ƃ��Č��o
  //
#if 0
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      if(eig0[i][j] >= thresh1) continue;
      if(eig0[i][j] >= thresh2 && eig1[i][j] >= thresh2){
        tmp[i][j] = 255;
      }
    }
  }
  //�����𖞂����_�̂����C���͂ɂ������𖞂����_�����݂�����̃G�b�W�����Ƃ���
  //��{�����͏�Ɠ��l

#elseif 0
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  int window = 3;
  bool exist;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      exist = false;
      if (!(is_edge(eig0[i][j], eig1[i][j], thresh1, thresh2)))
        continue;
      for (int ii = -window; ii < window; ++ii){
        for (int jj = -window; jj < window; ++jj){
          int i_ = std::min(std::max(0, ii + i), height);
          int j_ = std::min(std::max(0, jj + j), width);
          if (is_edge(eig0[i_][j_], eig1[i_][j_], thresh1, thresh2)){
            exist = true;
            break;
          }
        }
        if (exist) break;
      }
      if (exist)
        tmp[i][j] = 255;
    }
  }
#else
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      if(eig0[i][j] >= thresh1) continue;
      if(eig0[i][j] + eig1[i][j] >= thresh2 * 3){
        tmp[i][j] = 255;
      }
    }
  }
#endif
  // �R�[�i�[�_���̈�i�����̓_�̏W���j�̓_�ɑ΂��A�e�_�̑傫���ŗL�l��
  // �������ŗL�l�̍����ł����Ȃ��_���Ǐ��œK�_�Ƃ��ċ��߂�B�i���̍œK��������j
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      if(!tmp[i][j]) continue;
      int i0=i, j0=j, flag = 1;
      while(flag) {
         flag = 0;
         double min_dif = eig0[i0][j0] - eig1[i0][j0];
         for(int k = -5; k <= 5; k++){  //�O��T,���E�T�̗̈��T��
           for(int m = -5; m <= 5; m++){
             if(k == 0 && m == 0) continue;
	     if(i0+k < 0 || i0+k >= height) continue;
	     if(j0+m < 0 || j0+m >= width) continue;
             if(!tmp[i0+k][j0+m])continue;

	     if( (eig0[i0+k][j0+m] - eig1[i0+k][j0+m]) < min_dif){
	          min_dif = eig0[i0+k][j0+m] - eig1[i0+k][j0+m];
	          i0 = i0+k; j0 = j0+m;
	          flag = 1;
	      }
	      else tmp[i0+k][j0+m] = 0;
           }
	  }
      }
      pic[i0][j0] = 255;
    }
  }
  //��Ɨ̈�̉��
  freeInt2D(tmp, height, width);
}
