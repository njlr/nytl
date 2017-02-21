Concept: Vector<T>
- Represents a mathematical vector of a finite-dimensional vector space.
- T represents the field of the Vector. Therefore the values of the vector are part of a field.
	- This means it must provide the +,*,-,/ operators, a
	  0 and 1 value and respect the associative/commutative properties of a field.
	- There must be a valid nytl::FieldTraits specialization for the field over which
	  the vector is defined.
- Every vector can be multiplied by a value of its field and be added to another vector
  of the same space. This operations must be implemented using the * and + operators.
- DefaultConstructible, CopyConstructable, CopyAssignable, Destructable, EqualityComparable.
- Provides a const size() function that returns the number of elements it stores aka
  the dimension this vector has.
- Provides const and non-const versions of the [] operator. Returns a mutable reference
  for the non-const version and a const-reference/copy for the const version.
  The operator can be called with an integer value that can hold the size of the vector.
- Must provide a template<Size, typename> Rebind member typedef for other vector types.
  If a Rebind with the given properies is not possible it must result in a compile-time error.
- Vectors are interpreted as column vectors regarding matrix operations.
- The size (or dimension) of a vector shall never be 0.


- __NOTE__: in comparison to a mathematical vector, the operations for a vector and
  the underlaying field don't have to have closure over their operations. This means
  there is no strict seperation between vectors of int and double and it's ok
  that the int vector/scalar will be converted to a double vector/scalar if needed.
  Otherwise int vectors would not be possible since the set described by the int
  type is not a field. Its enough if the resulting vectors/scalars for such operations are
  similiar. But this means that one has to e.g. implement multiplication with double for
  int vectors.


For an implementation of the concept see nytl/vec.
Pseudo-codish example:

```cpp
template<...>
class Vector {
public:
	using Size = ...; // usually std::size_t. Must be convertible from/to int.
	using Value = ...; // the value type, mathematical field
	using Reference = ...; // usually Value&
	using ConstReference = ...; // usually const Value&

	// rebinds the vector implementation class
	template<Size Dim, typename T> using Rebind = ...;

	static constexpr Size dim = ...; // dimension the vector has. Might be a symbolic value.
	Vector create(Size size); // creates a vector with the given size

public:
	Vector();
	~Vector();
	Vector(const Vector<T>&);
	Vector& operator=(const Vector<T>&);

	Size size() const;

	Reference operator[](Size);
	ConstReference operator[](Size) const;

	// Always needed are the begin and end function, i.e. the vector must be iteratable
	// [c|r]begin, [c|r]end,

	// Optionally all iterator typedefs and functions as in stl containers
	// data, max_size, at, fill
};

// The first operations must return a vector with the same dimension and over a similiar field.
auto operator*(Value, Vector);
auto operator+(Vector, Vector);
auto operator-(Vector, Vector);
bool operator==(Vector, Vector);
bool operator!=(Vector, Vector);
```
