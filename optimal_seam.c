#include <stdio.h>
#include <stdlib.h> 
#include <time.h>

// ############# Start section to remove #############
void create_i(double *i) {
	i[0] = 5; i[1] = 8; i[2] = 12; i[3] = 3;
	i[4] = 4; i[5] = 2; i[6] = 3; i[7] = 9;
	i[8] = 7; i[9] = 3; i[10] = 4; i[11] = 2;
	i[12] = 4; i[13] = 5; i[14] = 7; i[15] = 8;
}

void create_vertical_M(double *M) {
   M[0] = 5; M[1] = 8; M[2] = 12; M[3] = 3;
   M[4] = 9; M[5] = 7; M[6] = 6; M[7] = 12;
   M[8] = 14; M[9] = 9; M[10] = 10; M[11] = 8;
   M[12] = 14; M[13] = 14; M[14] = 15; M[15] = 16;
}

void create_horizontal_M(double *M) {
	M[0] = 5; M[1] = 12; M[2] = 18; M[3] = 12;
	M[4] = 4; M[5] = 6; M[6] = 9; M[7] = 18;
	M[8] = 7; M[9] = 7; M[10] = 10; M[11] = 11;
	M[12] = 4; M[13] = 9; M[14] = 14; M[15] = 18;
}

struct res_step3 {
	double vertical_optimal_cost;
	double horizontal_optimal_cost;
	double *vertical_M;
	double *horizontal_M;
};

struct res_step3 step3(int row, int col, double *i) {
	struct res_step3 res;
	res.vertical_optimal_cost = rand() % 50;
	res.horizontal_optimal_cost = rand() % 50;
	res.vertical_M = (double *)malloc(16 * sizeof(double));
	res.horizontal_M = (double *)malloc(16 * sizeof(double));
	create_vertical_M(res.vertical_M);
	create_horizontal_M(res.horizontal_M);
	return res;
}
// ############# End section to remove #############

// Data structure to hold information in a T cell
struct cell_T {
	// optimal cost to get to this dimension
	double optimal_cost;
	// cell coming from: 0 for up (horizontal seam),
	// 1 for left (vertical seam)
	int direction;
	// corresponding matrix for this dimension,
	// obtained based on optimal cost
	double *i;
};

int *backtrack_horizontal_seam(int row, int col, double *M) {
	int min_pos = col-1, crr_pos, next_min_pos;
	double min_value = M[min_pos];
	int crr_row = 0, next_row;

	int k;
	for (k = 2; k <= row; ++k) {
		crr_pos = k*col-1;
		if (M[crr_pos] < min_value) {
			min_value = M[crr_pos];
			min_pos = crr_pos;
			crr_row = k-1;
		}
	}

	int *path = (int *)malloc(col * sizeof(int));
	path[col-1] = min_pos;
	int crr_left = col-2;
	while(crr_left >= 0) {
		next_min_pos = min_pos-1;
		min_value = M[next_min_pos];
		// Check not first row.
		if (crr_row > 0) {
			crr_pos = min_pos-col-1;
			if (M[crr_pos] < min_value) {
				min_value = M[crr_pos];
				next_min_pos = crr_pos;
				next_row = crr_row-1;
			}
		}
		// Check not last row.
		if (crr_row < row-1) {
			crr_pos = min_pos+col-1;
			if (M[crr_pos] < min_value) {
				min_value = M[crr_pos];
				next_min_pos = crr_pos;
				next_row = crr_row-1;
			}
		}

		min_pos = next_min_pos;
		crr_row = next_row;
		path[crr_left] = min_pos;
		crr_left--;
	}

	return path;
}

int *backtrack_vertical_seam(int row, int col, double *M) {
	int min_pos = (row-1)*col, crr_pos, next_min_pos;
	double min_value = M[min_pos];
	int crr_col = 0, next_col;

	int k;
	for (k = 1; k < col; ++k) {
		crr_pos = (row-1)*col + k;
		if (M[crr_pos] < min_value) {
			min_value = M[crr_pos];
			min_pos = crr_pos;
			crr_col = k;
		}
	}

	int *path = (int *)malloc(row * sizeof(int));
	path[row-1] = min_pos;
	int crr_left = row-2;
	while(crr_left >= 0) {
		next_min_pos = min_pos-col;
		min_value = M[next_min_pos];
		// Check not first column.
		if (crr_col > 0) {
			crr_pos = min_pos-col-1;
			if (M[crr_pos] < min_value) {
				min_value = M[crr_pos];
				next_min_pos = crr_pos;
				next_col = crr_col-1;
			}
		}
		// Check not last column.
		if (crr_col < col-1) {
			crr_pos = min_pos-col+1;
			if (M[crr_pos] < min_value) {
				min_value = M[crr_pos];
				next_min_pos = crr_pos;
				next_col = crr_col+1;
			}
		}

		min_pos = next_min_pos;
		crr_col = next_col;
		path[crr_left] = min_pos;
		crr_left--;
	}

	return path;
}

void set_T_cell(int T_index, int T_prev_index,
	double move_cost, int direction, int m_rows,
	int m_cols, double *M, struct cell_T *T) {
	// Calculate optimal cost.
	T[T_index].optimal_cost = T[T_prev_index].optimal_cost + move_cost;

	// Set direction to represent cell above.
	T[T_index].direction = direction;

	int *path;
	if (direction == 0) {
		path = backtrack_horizontal_seam(m_rows, m_cols, M);

		// ############# Start print section #############
		int j;
		for (j = 0; j < m_cols; ++j)
			printf("%d ", path[j]);
		printf("\n");
		// ############# End print section #############

		// Build new image matrix for this cell.
		T[T_index].i = (double *)malloc((m_rows-1) * m_cols * sizeof(double));
		int k, l;
		for (k = 0; k < m_cols; ++k) { // construct each column at a time
			int crr_row = 0;
			for (l = 0; l < m_rows; ++l)
				if (l*m_cols + k != path[k]) { // check for elem to remove
					T[T_index].i[crr_row*m_cols+k]
						= T[T_prev_index].i[l*m_cols + k];
					crr_row++;
				}
		}

		// ############# Start print section #############
		for (k = 0; k < m_rows-1; ++k) {
			for (l = 0; l < m_cols; ++l)
				printf("%lf ", T[T_index].i[k*m_cols+l]);
			printf("\n");
		}
		// ############# End print section #############
	} else {
		path = backtrack_vertical_seam(m_rows, m_cols, M);

		// ############# Start print section #############
		int j;
		for (j = 0; j < m_rows; ++j)
			printf("%d ", path[j]);
		printf("\n");
		// ############# End print section #############

		// Build new image matrix for this cell.
		T[T_index].i = (double *)malloc(m_rows * (m_cols-1) * sizeof(double));
		int k, l;
		for (k = 0; k < m_rows; ++k) { // construct each row at a time
			int crr_col = 0;
			for (l = 0; l < m_cols; ++l)
				if (k*m_cols + l != path[k]) { // check for elem to remove
					T[T_index].i[k*m_cols+crr_col]
						= T[T_prev_index].i[k*m_cols + l];
					crr_col++;
				}
		}

		// ############# Start print section #############
		for (k = 0; k < m_rows; ++k) {
			for (l = 0; l < m_cols-1; ++l)
				printf("%lf ", T[T_index].i[k*m_cols+l]);
			printf("\n");
		}
		// ############# End print section #############
	}
}

void calculate_T_cell(int row_pos, int col_pos, int row, int col,
	struct cell_T *T) {
	int T_index = row_pos*col+col_pos;
	T[T_index].optimal_cost = -1;
	struct res_step3 res_horizontal;
	res_horizontal.horizontal_optimal_cost = -1;
	struct res_step3 res_vertical;
	res_vertical.vertical_optimal_cost = -1;
	int T_above_index, i_above_rows, i_above_cols;
	int T_left_index, i_left_rows, i_left_cols;

	printf("<%d %d>\n", row_pos, col_pos);

	// Check if first row, otherwise can look above.
	if (row_pos > 0) {
		// This computation reduces number of rows.
		T_above_index = (row_pos-1)*col+col_pos;
		i_above_rows = row-row_pos+1;
		i_above_cols = col-col_pos;
		res_horizontal = step3(i_above_rows, i_above_cols,
			T[T_above_index].i);
	}
	// Check if first column, otherwise can look left.
	if (col_pos > 0) {
		// This computation reduces the number of columns.
		T_left_index = T_index-1;
		i_left_rows = row-row_pos;
		i_left_cols = col-col_pos+1;
		res_vertical = step3(i_left_rows, i_left_cols,
			T[T_left_index].i);
	}

	if (res_horizontal.horizontal_optimal_cost != -1
		&& res_vertical.vertical_optimal_cost == -1) {
		// First column.
		set_T_cell(T_index, T_above_index,
			res_horizontal.horizontal_optimal_cost, 0,
			i_above_rows, i_above_cols,
			res_horizontal.horizontal_M, T);
	} else if (res_horizontal.horizontal_optimal_cost == -1
		&& res_vertical.vertical_optimal_cost != -1) {
		// First row.
		set_T_cell(T_index, T_left_index,
			res_vertical.vertical_optimal_cost, 1,
			i_left_rows, i_left_cols,
			res_vertical.vertical_M, T);
	} else {
		if (res_horizontal.horizontal_optimal_cost
			< res_vertical.vertical_optimal_cost) {
			// Optimal to take value above.
			set_T_cell(T_index, T_above_index,
				res_horizontal.horizontal_optimal_cost, 0,
				i_above_rows, i_above_cols,
				res_horizontal.horizontal_M, T);
		} else {
			// Optimal to take value left.
			set_T_cell(T_index, T_left_index,
				res_vertical.vertical_optimal_cost, 1,
				i_left_rows, i_left_cols,
				res_vertical.vertical_M, T);
		}
	}
}

struct cell_T *optimal_seam(int row, int col, double *i) {
	struct cell_T *T = (struct cell_T *)malloc(row * col * sizeof(struct cell_T));
	T[0].optimal_cost = 0;
	T[0].i = i;

	int j, k;
	for (j = 1; j < col; ++j)
		calculate_T_cell(0, j, row, col, T);
	for (j = 1; j < row; ++j)
		for (k = 0; k < col; ++k)
			calculate_T_cell(j, k, row, col, T);
}

int main() {
	srand(time(NULL));
	
	double *i = (double *)malloc(16 * sizeof(double));
	create_i(i);
	optimal_seam(4, 4, i);

	return 0;
}