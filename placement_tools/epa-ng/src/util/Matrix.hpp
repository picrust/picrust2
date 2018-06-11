#pragma once

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief
 *
 * @file
 * @ingroup utils
 */

#include <stdexcept>
#include <vector>

// =================================================================================================
//     Matrix
// =================================================================================================

template <typename T>
class Matrix
{
public:

    // -------------------------------------------------------------
    //     Typedefs
    // -------------------------------------------------------------

    using value_type    = T;
    using iterator      = typename std::vector<T>::iterator;
    using const_iterator= typename std::vector<T>::const_iterator;

    // -------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------

    Matrix()
        : rows_(0)
        , cols_(0)
        , data_()
    {}

    Matrix (size_t rows, size_t cols)
        : rows_(rows)
        , cols_(cols)
        , data_(rows * cols)
    {}

    Matrix (size_t rows, size_t cols, T init)
        : rows_(rows)
        , cols_(cols)
        , data_(rows * cols, init)
    {}

    Matrix (size_t rows, size_t cols, std::initializer_list<T> const& init_list)
        : rows_(rows)
        , cols_(cols)
        , data_(rows * cols)
    {
        if (init_list.size() != size()) {
            throw std::length_error("__FUNCTION__: length_error");
        }

        size_t i = 0;
        for (T const& v : init_list) {
            data_[i] = v;
            ++i;
        }
    }

    ~Matrix() = default;

    Matrix(Matrix const&) = default;
    Matrix(Matrix&&)      = default;

    Matrix& operator= (Matrix const&) = default;
    Matrix& operator= (Matrix&&)      = default;

    void swap (Matrix& other)
    {
        using std::swap;
        swap(rows_, other.rows_);
        swap(cols_, other.cols_);
        swap(data_, other.data_);
    }

    // -------------------------------------------------------------
    //     Properties
    // -------------------------------------------------------------

    size_t rows() const
    {
        return rows_;
    }

    size_t cols() const
    {
        return cols_;
    }

    size_t size() const
    {
        return rows_ * cols_;
    }

    // -------------------------------------------------------------
    //     Element Access
    // -------------------------------------------------------------

    inline size_t coord(const size_t row, const size_t col) const
    {
        return row * cols_ + col;
    }

    T& at (const size_t row, const size_t col)
    {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("__FUNCTION__: out_of_range");
        }
        return data_[row * cols_ + col];
    }

    const T at (const size_t row, const size_t col) const
    {
        if (row >= rows_ || col >= cols_) {
            throw std::out_of_range("__FUNCTION__: out_of_range");
        }
        return data_[row * cols_ + col];
    }

    T& operator () (const size_t row, const size_t col)
    {
        return data_[row * cols_ + col];
    }

    inline const T operator () (const size_t row, const size_t col) const
    {
        return data_[row * cols_ + col];
    }

    const std::vector<T>& get_array() const
    {
        return data_;
    }

    // -------------------------------------------------------------
    //     Slicing
    // -------------------------------------------------------------

    std::vector<T> row( size_t index ) const
    {
        if( index >= rows_ ) {
            throw std::out_of_range("__FUNCTION__: out_of_range");
        }

        auto result = std::vector<T>( cols() );
        for( size_t i = 0; i < cols(); ++i ) {
            result[i] = operator()( index, i );
        }
        return result;
    }

    std::vector<T> col( size_t index ) const
    {
        if( index >= cols_ ) {
            throw std::out_of_range("__FUNCTION__: out_of_range");
        }

        auto result = std::vector<T>( rows() );
        for( size_t i = 0; i < rows(); ++i ) {
            result[i] = operator()( i, index );
        }
        return result;
    }

    // -------------------------------------------------------------
    //     Iterators
    // -------------------------------------------------------------

    iterator begin()
    {
        return data_.begin();
    }

    iterator end()
    {
        return data_.end();
    }

    const_iterator begin() const
    {
        return data_.begin();
    }

    const_iterator end() const
    {
        return data_.end();
    }

    const_iterator cbegin() const
    {
        return data_.cbegin();
    }

    const_iterator cend() const
    {
        return data_.cend();
    }

    // -------------------------------------------------------------
    //     Operators
    // -------------------------------------------------------------

    bool operator == (const Matrix<T>& rhs) const
    {
        return rows_ == rhs.rows_
            && cols_ == rhs.cols_
            && data_ == rhs.data_;
    }

    bool operator != (const Matrix<T>& rhs) const
    {
        return !(*this == rhs);
    }

    // -------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------

private:

    size_t         rows_;
    size_t         cols_;
    std::vector<T> data_;
};
