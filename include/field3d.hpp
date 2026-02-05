#pragma once
#include <vector>
#include <cassert>
#include <algorithm>
#include <cstring>

/**
 * Field3D - Flat 3D array wrapper for cache-efficient memory layout
 * 
 * Replaces nested std::vector<std::vector<std::vector<float>>> with
 * a single contiguous array, improving cache locality and reducing memory overhead.
 * 
 * Memory layout: Row-major order (r varies slowest, z fastest)
 * Index calculation: idx = i * NTH * NZ + j * NZ + k
 */
class Field3D {
public:
    // Default constructor - creates empty field
    Field3D() : NR_(0), NTH_(0), NZ_(0) {}
    
    // Constructor with dimensions
    Field3D(int nr, int nth, int nz) : NR_(nr), NTH_(nth), NZ_(nz) {
        data_.resize(nr * nth * nz, 0.0f);
    }
    
    // Constructor with dimensions and initial value
    Field3D(int nr, int nth, int nz, float init_value) : NR_(nr), NTH_(nth), NZ_(nz) {
        data_.resize(nr * nth * nz, init_value);
    }
    
    // Copy constructor
    Field3D(const Field3D& other) : NR_(other.NR_), NTH_(other.NTH_), NZ_(other.NZ_), data_(other.data_) {}
    
    // Move constructor
    Field3D(Field3D&& other) noexcept : NR_(other.NR_), NTH_(other.NTH_), NZ_(other.NZ_), data_(std::move(other.data_)) {
        other.NR_ = other.NTH_ = other.NZ_ = 0;
    }
    
    // Copy assignment
    Field3D& operator=(const Field3D& other) {
        if (this != &other) {
            NR_ = other.NR_;
            NTH_ = other.NTH_;
            NZ_ = other.NZ_;
            data_ = other.data_;
        }
        return *this;
    }
    
    // Move assignment
    Field3D& operator=(Field3D&& other) noexcept {
        if (this != &other) {
            NR_ = other.NR_;
            NTH_ = other.NTH_;
            NZ_ = other.NZ_;
            data_ = std::move(other.data_);
            other.NR_ = other.NTH_ = other.NZ_ = 0;
        }
        return *this;
    }
    
    // Resize the field
    void resize(int nr, int nth, int nz) {
        NR_ = nr;
        NTH_ = nth;
        NZ_ = nz;
        data_.resize(nr * nth * nz, 0.0f);
    }
    
    // Resize with initial value
    void resize(int nr, int nth, int nz, float init_value) {
        NR_ = nr;
        NTH_ = nth;
        NZ_ = nz;
        data_.resize(nr * nth * nz, init_value);
    }
    
    // Fill entire field with a value
    void fill(float value) {
        std::fill(data_.begin(), data_.end(), value);
    }
    
    // Assign from nested vector (for migration compatibility)
    void assign(int nr, int nth, int nz, const std::vector<std::vector<std::vector<float>>>& nested) {
        resize(nr, nth, nz);
        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nth; ++j) {
                for (int k = 0; k < nz; ++k) {
                    (*this)[i][j][k] = nested[i][j][k];
                }
            }
        }
    }
    
    // Size queries
    int size_r() const { return NR_; }
    int size_th() const { return NTH_; }
    int size_z() const { return NZ_; }
    size_t size() const { return data_.size(); }
    bool empty() const { return data_.empty(); }
    
    // Raw data access (for SIMD/GPU)
    float* data() { return data_.data(); }
    const float* data() const { return data_.data(); }
    
    // Direct index access (for performance-critical code)
    float& operator()(int i, int j, int k) {
        assert(i >= 0 && i < NR_ && j >= 0 && j < NTH_ && k >= 0 && k < NZ_);
        return data_[i * NTH_ * NZ_ + j * NZ_ + k];
    }
    
    const float& operator()(int i, int j, int k) const {
        assert(i >= 0 && i < NR_ && j >= 0 && j < NTH_ && k >= 0 && k < NZ_);
        return data_[i * NTH_ * NZ_ + j * NZ_ + k];
    }
    
    // Proxy classes for [i][j][k] syntax
private:
    // Proxy for [i] access - stores i, returns proxy for [j]
    class Slice2D {
        Field3D* parent_;
        int i_;
        
    public:
        Slice2D(Field3D* parent, int i) : parent_(parent), i_(i) {}
        
        // Proxy for [j] access - stores i and j, returns proxy for [k]
        class Slice1D {
            Field3D* parent_;
            int i_;
            int j_;
            
        public:
            Slice1D(Field3D* parent, int i, int j) : parent_(parent), i_(i), j_(j) {}
            
            // Proxy for [k] access - stores i, j, k, returns reference
            class Slice0D {
                Field3D* parent_;
                int i_;
                int j_;
                int k_;
                
            public:
                Slice0D(Field3D* parent, int i, int j, int k) : parent_(parent), i_(i), j_(j), k_(k) {}
                
                // Implicit conversion to float& for assignment
                operator float&() {
                    return (*parent_)(i_, j_, k_);
                }
                
                // Assignment operator
                float& operator=(float value) {
                    return (*parent_)(i_, j_, k_) = value;
                }
                
                // Addition assignment
                float& operator+=(float value) {
                    return (*parent_)(i_, j_, k_) += value;
                }
                
                // Addition assignment for double (converts to float)
                float& operator+=(double value) {
                    return (*parent_)(i_, j_, k_) += static_cast<float>(value);
                }
                
                // Subtraction assignment
                float& operator-=(float value) {
                    return (*parent_)(i_, j_, k_) -= value;
                }
                
                // Subtraction assignment for double (converts to float)
                float& operator-=(double value) {
                    return (*parent_)(i_, j_, k_) -= static_cast<float>(value);
                }
                
                // Multiplication assignment
                float& operator*=(float value) {
                    return (*parent_)(i_, j_, k_) *= value;
                }
                
                // Division assignment
                float& operator/=(float value) {
                    return (*parent_)(i_, j_, k_) /= value;
                }
                
                // Comparison operators
                bool operator==(float value) const {
                    return (*parent_)(i_, j_, k_) == value;
                }
                
                bool operator!=(float value) const {
                    return (*parent_)(i_, j_, k_) != value;
                }
                
                // Conversion to float for reading
                operator float() const {
                    return (*parent_)(i_, j_, k_);
                }
            };
            
            // Access [k]
            Slice0D operator[](int k) {
                return Slice0D(parent_, i_, j_, k);
            }
            
            // Const version for [k]
            class ConstSlice0D {
                const Field3D* parent_;
                int i_;
                int j_;
                int k_;
                
            public:
                ConstSlice0D(const Field3D* parent, int i, int j, int k) : parent_(parent), i_(i), j_(j), k_(k) {}
                
                operator float() const {
                    return (*parent_)(i_, j_, k_);
                }
                
                // Explicit conversion to double (no implicit conversion to avoid ambiguity)
                double to_double() const {
                    return static_cast<double>((*parent_)(i_, j_, k_));
                }
            };
            
            ConstSlice0D operator[](int k) const {
                return ConstSlice0D(parent_, i_, j_, k);
            }
        };
        
        // Access [j]
        Slice1D operator[](int j) {
            return Slice1D(parent_, i_, j);
        }
        
        // Const version for [j]
        class ConstSlice1D {
            const Field3D* parent_;
            int i_;
            int j_;
            
        public:
            ConstSlice1D(const Field3D* parent, int i, int j) : parent_(parent), i_(i), j_(j) {}
            
            class ConstSlice0D {
                const Field3D* parent_;
                int i_;
                int j_;
                int k_;
                
            public:
                ConstSlice0D(const Field3D* parent, int i, int j, int k) : parent_(parent), i_(i), j_(j), k_(k) {}
                
                operator float() const {
                    return (*parent_)(i_, j_, k_);
                }
                
                // Explicit conversion to double (no implicit conversion to avoid ambiguity)
                double to_double() const {
                    return static_cast<double>((*parent_)(i_, j_, k_));
                }
            };
            
            ConstSlice0D operator[](int k) const {
                return ConstSlice0D(parent_, i_, j_, k);
            }
        };
        
        ConstSlice1D operator[](int j) const {
            return ConstSlice1D(parent_, i_, j);
        }
    };
    
    // Const version of Slice2D
    class ConstSlice2D {
        const Field3D* parent_;
        int i_;
        
    public:
        ConstSlice2D(const Field3D* parent, int i) : parent_(parent), i_(i) {}
        
        class ConstSlice1D {
            const Field3D* parent_;
            int i_;
            int j_;
            
        public:
            ConstSlice1D(const Field3D* parent, int i, int j) : parent_(parent), i_(i), j_(j) {}
            
            class ConstSlice0D {
                const Field3D* parent_;
                int i_;
                int j_;
                int k_;
                
            public:
                ConstSlice0D(const Field3D* parent, int i, int j, int k) : parent_(parent), i_(i), j_(j), k_(k) {}
                
                operator float() const {
                    return (*parent_)(i_, j_, k_);
                }
                
                // Explicit conversion to double (no implicit conversion to avoid ambiguity)
                double to_double() const {
                    return static_cast<double>((*parent_)(i_, j_, k_));
                }
            };
            
            ConstSlice0D operator[](int k) const {
                return ConstSlice0D(parent_, i_, j_, k);
            }
        };
        
        ConstSlice1D operator[](int j) const {
            return ConstSlice1D(parent_, i_, j);
        }
    };
    
public:
    // Access [i][j][k]
    Slice2D operator[](int i) {
        return Slice2D(this, i);
    }
    
    // Const access [i][j][k]
    ConstSlice2D operator[](int i) const {
        return ConstSlice2D(this, i);
    }
    
private:
    int NR_, NTH_, NZ_;
    std::vector<float> data_;
};
