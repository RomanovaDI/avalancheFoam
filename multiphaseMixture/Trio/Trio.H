/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Trio

Description
    An ordered trio of two objects of type \<T\> with first() and second()
    elements.

SourceFiles
    TrioI.H

See also
    Foam::Tuple2 for storing two objects of dissimilar types.

\*---------------------------------------------------------------------------*/

#ifndef Trio_H
#define Trio_H

#include "FixedList.H"
#include "Istream.H"
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<typename T1, typename T2, typename T3>
using trio = std::tuple<T1, T2, T3>;

namespace Foam
{

// Forward Declarations
template<class T> class Trio;

// Common trio types
typedef Trio<label> labelTrio;      //!< A trio of labels
typedef Trio<word> wordTrio;        //!< A trio of words


/*---------------------------------------------------------------------------*\
                            Class Trio Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Trio
:
    public FixedList<T, 3>
{
public:

    // Constructors

        //- Default construct
        Trio() = default;

        //- Copy construct from components
        inline Trio(const T& f, const T& s, const T& t);

        //- Move construct from components
        inline Trio(T&& f, T&& s, T&& t);

        //- Copy construct from std::trio
        inline Trio(const trio<T,T,T>& vals);

        //- Move construct from std::trio
        inline Trio(trio<T,T,T>&& vals);

        //- Copy construct FixedList of two items
        inline Trio(const FixedList<T, 3>& list);

        //- Copy construct, optionally sorted with first less-than second
        //inline Trio(const T& f, const T& s, const T& t, const bool doSort);

        //- Copy construct, optionally sorted with first less-than second
        //inline Trio(const FixedList<T, 3>& list, const bool doSort);

        //- Construct from Istream
        inline explicit Trio(Istream& is);


    // Member Functions

    // Access

        //- Return first element
        using FixedList<T, 3>::first;

        //- Return middle element
        //using FixedList<T, 3>::middle;

        //- Return last element
        using FixedList<T, 3>::last;
        
		using FixedList<T, 3>::cdata;
		using FixedList<T, 3>::data;

		//- Return element of FixedList
        //inline T& operator[](const label i);

        //- Return element of constant FixedList
        //inline const T& operator[](const label i) const;

        //- Return second element, which is also the middle element
        inline const T& second() const noexcept;

        //- Return second element, which is also the middle element
        inline T& second() noexcept;

        //- Return previous element
        inline const T& previous(const T& a) const;

        //- Return next element
        inline const T& next(const T& a) const;

        //- Return third element, which is also the last element
        inline const T& third() const noexcept;

        //- Return third element, which is also the last element
        inline T& third() noexcept;


    // Queries

        //- True if first() is less-than second()
        //inline bool sorted() const;


    // Editing

        //- Flip the Trio in-place.
        //inline void flip();

        //- Sort so that first() is less-than second()
        //inline void sort();


    // Comparison

        //- Compare Trios
        //  \return
        //  -  0: different
        //  - +1: identical values and order used
        //  - -1: identical values, but in reversed order
        static inline int compare(const Trio<T>& a, const Trio<T>& b);


    // Hashing

        //- Hashing functor for Trio, hashes lower value first
        //  Regular hasher inherited from FixedList
        /*struct hasher
        {
            unsigned operator()(const Trio<T>& obj, unsigned seed=0) const
            {
                Foam::Hash<T> op;
                if (obj.second() < obj.first())
                {
                    return op(obj.first(), op(obj.second(), seed));
                }
                else
                {
                    return op(obj.second(), op(obj.first(), seed));
                }
            }
        };*/
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Trio is contiguous if the type is contiguous
template<class T>
struct is_contiguous<Trio<T>> : is_contiguous<T> {};

//- Check for Trio of labels
template<class T>
struct is_contiguous_label<Trio<T>> : is_contiguous_label<T> {};

//- Check for Trio of scalars
template<class T>
struct is_contiguous_scalar<Trio<T>> : is_contiguous_scalar<T> {};

//- Hashing for Trio of data
template<class T>
struct Hash<Trio<T>> : Trio<T>::hasher {};


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return reverse of a Trio
/*template<class T>
Trio<T> reverse(const Trio<T>& p)
{
    return Trio<T>(p.third(), p.second(), p.first());
}*/


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

template<class T>
bool operator==(const Trio<T>& a, const Trio<T>& b)
{
    return (a.first() == b.first() && a.second() == b.second() && a.third() == b.third());
}


template<class T>
bool operator!=(const Trio<T>& a, const Trio<T>& b)
{
    return !(a == b);
}


template<class T>
bool operator<(const Trio<T>& a, const Trio<T>& b)
{
    return
    (
        a.first() < b.first()
     || (!(b.first() < a.first()) && a.second() < b.second())
     || (!(b.first() < a.first()) && !(b.second() < a.second()) && a.third() < b.third())
    );
}


template<class T>
bool operator<=(const Trio<T>& a, const Trio<T>& b)
{
    return !(b < a);
}


template<class T>
bool operator>(const Trio<T>& a, const Trio<T>& b)
{
    return (b < a);
}


template<class T>
bool operator>=(const Trio<T>& a, const Trio<T>& b)
{
    return !(a < b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TrioI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
