// signature(obs="numeric", n_contigs="integer", n_states="integer", result="numeric")

#define INDEX(n_rows,row_i,col_i) col_i * n_rows + row_i
#define MIN(a,b) a<b ? a : b

int nn_contigs = *n_contigs;
int nn_states = *n_states;

for (unsigned int i1 = 0; i1 < nn_contigs; ++i1)
  for (unsigned int i2 = 0; i2 < nn_contigs; ++i2) {
    double a_i1=0, a_i2=0, a_both=0;
    for (int j=0; j<nn_states; j++) {
      a_i1 += obs[INDEX(nn_contigs,i1,j)] * obs[INDEX(nn_contigs,i1,j)];
      a_i2 += obs[INDEX(nn_contigs,i2,j)] * obs[INDEX(nn_contigs,i2,j)];
      a_both += obs[INDEX(nn_contigs,i1,j)] * obs[INDEX(nn_contigs,i2,j)];
    }
    double div = sqrt(a_i2) * sqrt(a_i1);
    if (div == 0)
      result[INDEX(nn_contigs,i1,i2)] = -1;
    else
      result[INDEX(nn_contigs,i1,i2)] = (a_both / div);
  }
