#include "auxFunctions.h"
#include <cctype> // This header declares a set of functions to classify and transform individual characters, like toupper
#include <cwctype> // Same thing for wide characters

/*****************************************************
    Useful math/geometric operations
*****************************************************/
mat::fixed<3, 3> rotatMatrix(const double rotatX, const double rotatY, const double rotatZ)
{
	/* Rotation matrix is RotatMat = RotatX * RotatY * RotatZ, i.e. a rotation around the Z axis,
	  followed by a rotation around the new y axis, and a rotation around the new x axis. Each rotation matrix is given by:
	
	  rotatX = { 	{ 1 ,        0        ,         0        },
			     			{ 0 , cos(rotation(3)) , -sin(rotation(3)) },
			     			{ 0 , sin(rotation(3)) ,  cos(rotation(3)) }
			   			};

	  rotatY = { 	{ cos(rotation(4))  , 0 ,  sin(rotation(4)) },
		         		{        0         , 1 ,         0        },
			     			{ -sin(rotation(4)) , 0 , cos(rotation(4)) }
			   			};

	  rotatZ = {	{ cos(rotation(5)) , -sin(rotation(5)) , 0 },			     
				 				{ sin(rotation(5)) ,  cos(rotation(5)) , 0 },
			     			{        0        ,         0        , 1 },
			   			};

	  And the resulting matrix is the one below
	*/
	mat::fixed<3, 3> rotatMatrix = { 
									{                          cos(rotatY) * cos(rotatZ)                               ,                          -cos(rotatY) * sin(rotatZ)                               ,            sin(rotatY)          },
									{ cos(rotatX) * sin(rotatZ) + sin(rotatX) * sin(rotatY) * cos(rotatZ)  ,  cos(rotatX) * cos(rotatZ) - sin(rotatX) * sin(rotatY) * sin(rotatZ)  ,  -sin(rotatX) * cos(rotatY) },
									{ sin(rotatX) * sin(rotatZ) - cos(rotatX) * sin(rotatY) * cos(rotatZ)  ,  sin(rotatX) * cos(rotatZ) + cos(rotatX) * sin(rotatY) * sin(rotatZ)  ,   cos(rotatX) * cos(rotatY) }
								   };	

	return rotatMatrix;
}

mat::fixed<3, 3> rotatMatrix_deg(const double rotatX, const double rotatY, const double rotatZ)
{
	return rotatMatrix(deg2rad(rotatX), deg2rad(rotatY), deg2rad(rotatZ));
}


mat::fixed<3, 3> rotatMatrix(const vec::fixed<3> &rotation)
{
	return rotatMatrix(rotation(0),rotation(1), rotation(2));
}

mat::fixed<3, 3> rotatMatrix_deg(const vec::fixed<3> &rotation)
{
	return rotatMatrix_deg(rotation(0), rotation(1), rotation(2));
}

double deg2rad(const double degree)
{
	return degree * datum::pi / 180;
}

double minimum(const double x, const double y)
{
	return (x > y) ? y : x;
}

/*****************************************************
    String handling
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
bool thereIsCommentInString(const std::string& str)
{
	return (str.find("%") != std::string::npos);
}

// Verify whether a string has content, i.e. if:
// 1) it is not empty
// 2) it is not just white spaces or tabs
// 3) it does not start with a comment mark ('%')
bool hasContent(const std::string &str)
{
	// Empty strings have no content
	if (str.empty())
		return false;

	// The string has content only if it has at least one character that is neither a white-space nor a tab (\t)
	return ((str.find_first_not_of(" \t") != std::string::npos) && (str.at(0) != '%'));
}

void removeComments(std::string &str)
{
	str = str.substr(0, str.find("%", 0));
}

// Tokenize a string using a given delimiter.
// Return a std::vector with the resulting strings at the different elements
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim)
{
	std::vector<std::string> tokens;
	std::string aux = str;

	while (hasContent(aux))
	{
		if (aux.find_first_not_of(delim) == std::string::npos) // Check if there is something besides delimiters in the line
		{
			break; // Then we break the while loop and return tokens as an empty std::vector
		}

		aux = aux.substr(aux.find_first_not_of(delim)); // Get the part of aux starting at the first character that is not a delimiter
		tokens.push_back(aux.substr(0, aux.find_first_of(delim))); // Take content before next delimiter and add to tokens

		if (aux.find_first_of(delim, 0) != std::string::npos) // If there is another delimiter in the string...
		{
			aux = aux.substr(aux.find_first_of(delim, 0)); // ... keep in aux everything after the second delimiter, including the delimiter
		}
		else // Otherwise, we have already included the last element in tokens, so we can end this loop.
		{
			break;
		}
	}

	return tokens;
}


// Convert lowercase letter to uppercase and compare if they are equal
inline bool caseInsCharCompareN(char a, char b) {
	return(std::toupper(a) == std::toupper(b));
}

// Same thing for wchar
inline bool caseInsCharCompareW(wchar_t a, wchar_t b) {
	return(std::towupper(a) == std::towupper(b));
}

// Compare each character of the strings to see if they match
bool caseInsCompare(const std::string& s1, const std::string& s2) {
	return((s1.size() == s2.size()) &&
		std::equal(s1.begin(), s1.end(), s2.begin(), caseInsCharCompareN));
}

// Same thing for wstring
bool caseInsCompare(const std::wstring& s1, const std::wstring& s2) {
	return((s1.size() == s2.size()) &&
		std::equal(s1.begin(), s1.end(), s2.begin(), caseInsCharCompareW));
}


// Get folder path from a complete file path
std::string getFileFolder(const std::string& path)
{
	// Check if input string is empty
	if (path.empty())
	{
		throw std::runtime_error("Empty string passed to getFileFolder().");
	}

	std::vector<std::string> str_tokenized = stringTokenize(path, filesep);

	// If there is only one element in str_tokenized, it means that only the file name
	// was provided as an argument. Hence, the fileFolder is the current directory, "."
	if (str_tokenized.size() == 1)
	{
		return ".";
	}

	// Otherwise, return the full path until the last delimiter
	else
	{
		size_t found;
		found = path.find_last_of(filesep);
		return path.substr(0, found);
	}
}


// Get file name, without extension, from a complete file path
std::string getFileName(const std::string& path)
{
	// Check if input string is empty
	if (path.empty())
	{
		throw std::runtime_error("Empty string passed to getFileName().");
	}

	// Tokenize the passed string based on the file separator. The file name is the last part (works even if only the file name is passed as argument)
	std::vector<std::string> str_tokenized = stringTokenize(path, filesep);
	std::string flNm = str_tokenized.back();

	// We need to exclude the file extension. We get everything until the last dot and say that 
	if (flNm.find_last_not_of(".") != std::string::npos)
	{
		flNm = flNm.substr(0, flNm.find_last_of("."));
	}

	return flNm;
}


/*****************************************************
	Spline implementation
	Downloaded from https://kluge.in-chemnitz.de/opensource/spline/
*****************************************************/

/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */

namespace tk
{

	// band_matrix implementation
	// -------------------------

	band_matrix::band_matrix(int dim, int n_u, int n_l)
	{
		resize(dim, n_u, n_l);
	}
	void band_matrix::resize(int dim, int n_u, int n_l)
	{
		assert(dim > 0);
		assert(n_u >= 0);
		assert(n_l >= 0);
		m_upper.resize(n_u + 1);
		m_lower.resize(n_l + 1);
		for (size_t i = 0; i < m_upper.size(); i++) {
			m_upper[i].resize(dim);
		}
		for (size_t i = 0; i < m_lower.size(); i++) {
			m_lower[i].resize(dim);
		}
	}
	int band_matrix::dim() const
	{
		if (m_upper.size() > 0) {
			return m_upper[0].size();
		}
		else {
			return 0;
		}
	}

	int band_matrix::num_upper() const
	{
		return m_upper.size() - 1;
	}
	int band_matrix::num_lower() const
	{
		return m_lower.size() - 1;
	}


	// defines the new operator (), so that we can access the elements
	// by A(i,j), index going from i=0,...,dim()-1
	double & band_matrix::operator () (int i, int j)
	{
		int k = j - i;       // what band is the entry
		assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
		assert((-num_lower() <= k) && (k <= num_upper()));
		// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
		if (k >= 0)   return m_upper[k][i];
		else	    return m_lower[-k][i];
	}
	double band_matrix::operator () (int i, int j) const
	{
		int k = j - i;       // what band is the entry
		assert((i >= 0) && (i < dim()) && (j >= 0) && (j < dim()));
		assert((-num_lower() <= k) && (k <= num_upper()));
		// k=0 -> diogonal, k<0 lower left part, k>0 upper right part
		if (k >= 0)   return m_upper[k][i];
		else	    return m_lower[-k][i];
	}
	// second diag (used in LU decomposition), saved in m_lower
	double band_matrix::saved_diag(int i) const
	{
		assert((i >= 0) && (i < dim()));
		return m_lower[0][i];
	}
	double & band_matrix::saved_diag(int i)
	{
		assert((i >= 0) && (i < dim()));
		return m_lower[0][i];
	}

	// LR-Decomposition of a band matrix
	void band_matrix::lu_decompose()
	{
		int  i_max, j_max;
		int  j_min;
		double x;

		// preconditioning
		// normalize column i so that a_ii=1
		for (int i = 0; i < this->dim(); i++) {
			assert(this->operator()(i, i) != 0.0);
			this->saved_diag(i) = 1.0 / this->operator()(i, i);
			j_min = std::max(0, i - this->num_lower());
			j_max = std::min(this->dim() - 1, i + this->num_upper());
			for (int j = j_min; j <= j_max; j++) {
				this->operator()(i, j) *= this->saved_diag(i);
			}
			this->operator()(i, i) = 1.0;          // prevents rounding errors
		}

		// Gauss LR-Decomposition
		for (int k = 0; k < this->dim(); k++) {
			i_max = std::min(this->dim() - 1, k + this->num_lower());  // num_lower not a mistake!
			for (int i = k + 1; i <= i_max; i++) {
				assert(this->operator()(k, k) != 0.0);
				x = -this->operator()(i, k) / this->operator()(k, k);
				this->operator()(i, k) = -x;                         // assembly part of L
				j_max = std::min(this->dim() - 1, k + this->num_upper());
				for (int j = k + 1; j <= j_max; j++) {
					// assembly part of R
					this->operator()(i, j) = this->operator()(i, j) + x * this->operator()(k, j);
				}
			}
		}
	}
	// solves Ly=b
	std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const
	{
		assert(this->dim() == (int)b.size());
		std::vector<double> x(this->dim());
		int j_start;
		double sum;
		for (int i = 0; i < this->dim(); i++) {
			sum = 0;
			j_start = std::max(0, i - this->num_lower());
			for (int j = j_start; j < i; j++) sum += this->operator()(i, j)*x[j];
			x[i] = (b[i] * this->saved_diag(i)) - sum;
		}
		return x;
	}
	// solves Rx=y
	std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const
	{
		assert(this->dim() == (int)b.size());
		std::vector<double> x(this->dim());
		int j_stop;
		double sum;
		for (int i = this->dim() - 1; i >= 0; i--) {
			sum = 0;
			j_stop = std::min(this->dim() - 1, i + this->num_upper());
			for (int j = i + 1; j <= j_stop; j++) sum += this->operator()(i, j)*x[j];
			x[i] = (b[i] - sum) / this->operator()(i, i);
		}
		return x;
	}

	std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
		bool is_lu_decomposed)
	{
		assert(this->dim() == (int)b.size());
		std::vector<double>  x, y;
		if (is_lu_decomposed == false) {
			this->lu_decompose();
		}
		y = this->l_solve(b);
		x = this->r_solve(y);
		return x;
	}




	// spline implementation
	// -----------------------

	void spline::set_boundary(spline::bd_type left, double left_value,
		spline::bd_type right, double right_value,
		bool force_linear_extrapolation)
	{
		assert(m_x.size() == 0);          // set_points() must not have happened yet
		m_left = left;
		m_right = right;
		m_left_value = left_value;
		m_right_value = right_value;
		m_force_linear_extrapolation = force_linear_extrapolation;
	}


	void spline::set_points(const std::vector<double>& x,
		const std::vector<double>& y, bool cubic_spline)
	{
		assert(x.size() == y.size());
		assert(x.size() > 2);
		m_x = x;
		m_y = y;
		int   n = x.size();
		// TODO: maybe sort x and y, rather than returning an error
		for (int i = 0; i < n - 1; i++) {
			assert(m_x[i] < m_x[i + 1]);
		}

		if (cubic_spline == true) { // cubic spline interpolation
			// setting up the matrix and right hand side of the equation system
			// for the parameters b[]
			band_matrix A(n, 1, 1);
			std::vector<double>  rhs(n);
			for (int i = 1; i < n - 1; i++) {
				A(i, i - 1) = 1.0 / 3.0*(x[i] - x[i - 1]);
				A(i, i) = 2.0 / 3.0*(x[i + 1] - x[i - 1]);
				A(i, i + 1) = 1.0 / 3.0*(x[i + 1] - x[i]);
				rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
			}
			// boundary conditions
			if (m_left == spline::second_deriv) {
				// 2*b[0] = f''
				A(0, 0) = 2.0;
				A(0, 1) = 0.0;
				rhs[0] = m_left_value;
			}
			else if (m_left == spline::first_deriv) {
				// c[0] = f', needs to be re-expressed in terms of b:
				// (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
				A(0, 0) = 2.0*(x[1] - x[0]);
				A(0, 1) = 1.0*(x[1] - x[0]);
				rhs[0] = 3.0*((y[1] - y[0]) / (x[1] - x[0]) - m_left_value);
			}
			else {
				assert(false);
			}
			if (m_right == spline::second_deriv) {
				// 2*b[n-1] = f''
				A(n - 1, n - 1) = 2.0;
				A(n - 1, n - 2) = 0.0;
				rhs[n - 1] = m_right_value;
			}
			else if (m_right == spline::first_deriv) {
				// c[n-1] = f', needs to be re-expressed in terms of b:
				// (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
				// = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
				A(n - 1, n - 1) = 2.0*(x[n - 1] - x[n - 2]);
				A(n - 1, n - 2) = 1.0*(x[n - 1] - x[n - 2]);
				rhs[n - 1] = 3.0*(m_right_value - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
			}
			else {
				assert(false);
			}

			// solve the equation system to obtain the parameters b[]
			m_b = A.lu_solve(rhs);

			// calculate parameters a[] and c[] based on b[]
			m_a.resize(n);
			m_c.resize(n);
			for (int i = 0; i < n - 1; i++) {
				m_a[i] = 1.0 / 3.0*(m_b[i + 1] - m_b[i]) / (x[i + 1] - x[i]);
				m_c[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
					- 1.0 / 3.0*(2.0*m_b[i] + m_b[i + 1])*(x[i + 1] - x[i]);
			}
		}
		else { // linear interpolation
			m_a.resize(n);
			m_b.resize(n);
			m_c.resize(n);
			for (int i = 0; i < n - 1; i++) {
				m_a[i] = 0.0;
				m_b[i] = 0.0;
				m_c[i] = (m_y[i + 1] - m_y[i]) / (m_x[i + 1] - m_x[i]);
			}
		}

		// for left extrapolation coefficients
		m_b0 = (m_force_linear_extrapolation == false) ? m_b[0] : 0.0;
		m_c0 = m_c[0];

		// for the right extrapolation coefficients
		// f_{n-1}(x) = b*(x-x_{n-1})^2 + c*(x-x_{n-1}) + y_{n-1}
		double h = x[n - 1] - x[n - 2];
		// m_b[n-1] is determined by the boundary condition
		m_a[n - 1] = 0.0;
		m_c[n - 1] = 3.0*m_a[n - 2] * h*h + 2.0*m_b[n - 2] * h + m_c[n - 2];   // = f'_{n-2}(x_{n-1})
		if (m_force_linear_extrapolation == true)
			m_b[n - 1] = 0.0;
	}

	double spline::operator() (double x) const
	{
		size_t n = m_x.size();
		// find the closest point m_x[idx] < x, idx=0 even if x<m_x[0]
		std::vector<double>::const_iterator it;
		it = std::lower_bound(m_x.begin(), m_x.end(), x);
		int idx = std::max(int(it - m_x.begin()) - 1, 0);

		double h = x - m_x[idx];
		double interpol;
		if (x < m_x[0]) {
			// extrapolation to the left
			interpol = (m_b0*h + m_c0)*h + m_y[0];
		}
		else if (x > m_x[n - 1]) {
			// extrapolation to the right
			interpol = (m_b[n - 1] * h + m_c[n - 1])*h + m_y[n - 1];
		}
		else {
			// interpolation
			interpol = ((m_a[idx] * h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
		}
		return interpol;
	}

}
