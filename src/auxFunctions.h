#pragma once
#include <armadillo>
#include <string>

// These four below are needed for the spline implementation
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>


// Define file separator for current platform
const std::string filesep =
#ifdef _WIN32
	"\\";
#else
	"/";
#endif

typedef std::vector<std::complex<double>> cx_stdvec;


/*****************************************************
    Useful math/geometric operations
*****************************************************/
arma::mat::fixed<3, 3> rotatMatrix(const double rotatX, const double rotatY, const double rotatZ);
arma::mat::fixed<3, 3> rotatMatrix(const arma::vec::fixed<3> &rotation);
arma::mat::fixed<3, 3> rotatMatrix_deg(const double rotatX, const double rotatY, const double rotatZ);
arma::mat::fixed<3, 3> rotatMatrix_deg(const arma::vec::fixed<3> &rotation);
arma::mat::fixed<3, 3> rotatMatrix_extrinsic(const double rotatX, const double rotatY, const double rotatZ);
arma::mat::fixed<3, 3> rotatMatrix_extrinsic(const arma::vec::fixed<3> &rotation);
arma::mat::fixed<3, 3> smallRotatMatrix(const double rotatX, const double rotatY, const double rotatZ);
arma::mat::fixed<3, 3> smallRotatMatrix(const arma::vec::fixed<3> &rotation);
double deg2rad(const double degree);
double minimum(const double x, const double y);
bool almostEqual(const double x, const double y, double eps);


/*****************************************************
    String handling
*****************************************************/

// Verify whether a string contains a comment, marked by a '%'
// TODO: define comment character in the beginning of this file
bool thereIsCommentInString(const std::string& str);

// Verify whether a string has content, i.e. if it is not empty, it is not just
// white spaces or tabs, and it does not start with a comment mark ('%')
bool hasContent(const std::string &str);

// Remove comments from string, marked by a '%'
void removeComments(std::string &str);

// Tokenize a string using a given delimiter
std::vector<std::string> stringTokenize(const std::string &str, const std::string &delim = " \t");

// Case-insensitive string comparison
// Found at https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch04s14.html
// Convert lowercase letter to uppercase and compare if they are equal
inline bool caseInsCharCompareN(char a, char b);

// Same thing for wchar
inline bool caseInsCharCompareW(wchar_t a, wchar_t b);

// Compare each character of the strings to see if they match
bool caseInsCompare(const std::string& s1, const std::string& s2);

// Same thing for wstring
bool caseInsCompare(const std::wstring& s1, const std::wstring& s2);

// Get folder path from a complete file path
std::string getFileFolder(const std::string& path);

// Get file name, without extension, from a complete file path
std::string getFileName(const std::string& path);

// string2num: used to convert from string to a numerical type (double, float, int...)
// Returns True if the conversion is succesful and False if it is not
// Found at http://www.learncpp.com/cpp-tutorial/17-2-stdstring-construction-and-destruction/
template <typename T>
inline bool string2num(const std::string& sString, T &tX)
{
	std::istringstream iStream(sString);
	return (iStream >> tX) ? true : false; // extract value into tX, return success or not
}

// string2num: used to convert from string to a numerical type (double, float, int...)
// Adapted from the one found at http://www.learncpp.com/cpp-tutorial/17-2-stdstring-construction-and-destruction/
//
// Examples:
// --> double a = string2num<double>("15");
// leads to a == 15
//
// --> double a = estring2num<double>("15 32");
// --> double a = string2num<double>("15 a");
// --> double a = string2num<double>("15a32");
// --> double a = string2num<double>("15a");
// all lead to a == 15
//
// --> double a = string2num<double>("a");
// throws an exception
template <typename T>
inline T string2num(const std::string& string)
{
	std::istringstream iStream(string);
	T tX;
	if (!(iStream >> tX))
	{
		throw std::runtime_error( "Conversion from string failed. Tried to convert the string'" + string + "' to a number.");
	}
	return tX;
}


// FFT and IFFT functions
cx_stdvec mkl_ifft(cx_stdvec &in);
std::vector<double> mkl_ifft_real(cx_stdvec &in);
arma::mat mkl_ifft_real(arma::cx_mat &in);


/*****************************************************
	Spline implementation
	It is almost equal to the original file downloaded from https://kluge.in-chemnitz.de/opensource/spline/
	but it has some minor modifications
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

	// band matrix solver
	class band_matrix
	{
	private:
		std::vector< std::vector<double> > m_upper;  // upper band
		std::vector< std::vector<double> > m_lower;  // lower band
	public:
		band_matrix() {};                             // constructor
		band_matrix(int dim, int n_u, int n_l);       // constructor
		~band_matrix() {};                            // destructor
		void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
		int dim() const;                             // matrix dimension
		int num_upper() const;
		int num_lower() const;
		// access operator
		double & operator () (int i, int j);            // write
		double   operator () (int i, int j) const;      // read
		// we can store an additional diogonal (in m_lower)
		double& saved_diag(int i);
		double  saved_diag(int i) const;
		void lu_decompose();
		std::vector<double> r_solve(const std::vector<double>& b) const;
		std::vector<double> l_solve(const std::vector<double>& b) const;
		std::vector<double> lu_solve(const std::vector<double>& b,
			bool is_lu_decomposed = false);

	};


	// spline interpolation
	class spline
	{
	public:
		enum bd_type {
			first_deriv = 1,
			second_deriv = 2
		};

	private:
		std::vector<double> m_x, m_y;            // x,y coordinates of points
		// interpolation parameters
		// f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
		std::vector<double> m_a, m_b, m_c;        // spline coefficients
		double  m_b0, m_c0;                     // for left extrapol
		bd_type m_left, m_right;
		double  m_left_value, m_right_value;
		bool    m_force_linear_extrapolation;

	public:
		// set default boundary condition to be zero curvature at both ends
		spline() : m_left(second_deriv), m_right(second_deriv),
			m_left_value(0.0), m_right_value(0.0),
			m_force_linear_extrapolation(false)
		{
			;
		}

		// optional, but if called it has to come be before set_points()
		void set_boundary(bd_type left, double left_value,
			bd_type right, double right_value,
			bool force_linear_extrapolation = false);
		void set_points(const std::vector<double>& x,
			const std::vector<double>& y, bool cubic_spline = true);
		double operator() (double x) const;
	};
}
