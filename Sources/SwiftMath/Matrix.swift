//
//  Matrix.swift
//  MySwift
//
//  Created by Jacob W Esselstyn on 1/19/25.
//
import Foundation

public struct Matrix<T: BinaryFloatingPoint>: Hashable {
	
	public var rows: Int { values.count }
	public var columns: Int { values.first?.count ?? 0 }
	public private(set) var values: [[T]]
	
	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
	public init(_ values: [[T]]) {
		self.values = values
	}
	public init(repeating value: T, rows: Int, columns: Int) {
		self.values = Array(repeating: Array(repeating: value, count: columns), count: rows)
	}
	
	/// Check if rows and columns of each matrix match.
	public static func areDimensionsEqual(lhs: Self, rhs: Self) -> Bool {
		lhs.rows == rhs.rows && lhs.columns == rhs.columns
	}
	/// Perform Matrix Multiplication
	public static func * (lhs: Self, rhs: Self) throws -> Self {
		guard lhs.rows == rhs.columns && lhs.columns == rhs.rows else { throw MatrixError.nonMatchingDimensions }
		var result: [[T]] = Array<Array<T>>(repeating: Array<T>(repeating: T.init(0.0), count: rhs.columns), count: lhs.rows)
		for i in 0..<lhs.rows {
			for j in 0..<rhs.columns {
				for k in 0..<lhs.columns {
					result[i][j] += lhs.values[i][k] * rhs.values[k][j]
				}
			}
		}
		return Matrix(result)
	}
	
	/// Elementwise addition
	public static func + (lhs: Self, rhs: Self) throws -> Self {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in result.values.indices {
			for j in result.values[i].indices {
				result.values[i][j] = lhs.values[i][j] + rhs.values[i][j]
			}
		}
		return result
	}
	/// Elementwise subtraction
	public static func - (lhs: Self, rhs: Self) throws -> Self {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in result.values.indices {
			for j in result.values[i].indices {
				result.values[i][j] = lhs.values[i][j] - rhs.values[i][j]
			}
		}
		return result
	}
	
	/// Elementwise multiplication
	public static func * (lhs: Self, rhs: T) -> Self {
		var result = lhs
		for i in result.values.indices {
			for j in result.values[i].indices {
				result.values[i][j] = result.values[i][j]*rhs
			}
		}
		return result
	}
	/// Elementwise multiplication
	public static func * (lhs: T, rhs: Self) -> Self {
		var result = rhs
		for i in result.values.indices {
			for j in result.values[i].indices {
				result.values[i][j] = result.values[i][j]*lhs
			}
		}
		return result
	}
	
	/// Elementwise multiplication of two matricies
	public mutating func multiplyElementwise(by other: Self) {
		for i in values.indices {
			for j in values[i].indices {
				values[i][j] = values[i][j]*other.values[i][j]
			}
		}
	}
	/// Elementwise multiplication of two matricies
	public func multiplyingElementwise(by other: Self) -> Self {
		var newMatrix = self
		newMatrix.multiplyElementwise(by: other)
		return newMatrix
	}
	
	/// Returns the transpose of the matrix (swapping rows for columns)
	public var transpose: Self {
		guard rows > 0 && columns > 0 else { return self }
		var newMatrix: [[T]] = Array<Array<T>>(repeating: Array<T>(repeating: T.init(0.0), count: rows), count: columns)
		for i in self.values.indices {
			for j in self.values[i].indices {
				newMatrix[j][i] = self.values[i][j]
			}
		}
		return Matrix(newMatrix)
	}
	
	/// Returns `true` if matrix is symmetric
	public var isSymmetric: Bool { self == self.transpose }
	
	/// Check if the matrix is square.
	private var isSquare: Bool { rows == columns }
	
	/// Determinant of the matrix (only defined for square matrices).
	public var determinant: T? {
		guard isSquare else { return nil }
		return Matrix.calculateDeterminant(for: values)
		//return rowEcholonForm?.mainDiagonal?.reduce(1, *) ?? 0
	}
	
	/// Recursive helper function to calculate the determinant.
	private static func calculateDeterminant(for matrix: [[T]]) -> T {
		let n = matrix.count
	
		if n == 1 { // Base case for 1x1 matrix
			return matrix[0][0]
		}
		if n == 2 { // Base case for 2x2 matrix
			return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
		}
	
		var det: T = 0
		// Calculate determinant using row reduction (Gaussian elimination) - O(n^3)
		// TO-DO
		// Calculate determinant using cofactor expansion - O(n!)
		for i in 0..<n {
			// Create the submatrix by excluding the first row and the ith column
			var subMatrix = [[T]]()
			for j in 1..<matrix.count {
				var row = [T]()
				for k in 0..<matrix.count {
					if k != i {
						row.append(matrix[j][k])
					}
				}
				subMatrix.append(row)
			}
			// Recursive determinant calculation with cofactor expansion
			let sign: T = (i % 2 == 0) ? 1 : -1
			det += sign * matrix[0][i] * calculateDeterminant(for: subMatrix)
		}
		return det
	}
	
	/// Fetches or updates a value in the matrix. Will crash if out-of-range, and does nothing when setting a value if out-of-range.
	public subscript(i: Int, j: Int) -> T {
		get {
			if i >= rows { fatalError("row index out of range of Matrix") }
			if j >= columns { fatalError("column index out of range of Matrix") }
			return values[i][j]
		}
		set(newValue) {
			if i >= rows && j >= columns { return }
			values[i][j] = newValue
		}
	}
	/// Fetches or updates a value in the matrix. Returns `nil` if out-of-range, and does nothing when setting a value if out-of-range or if new value is `nil`.
	public subscript(i: Int, j: Int) -> T? {
		get {
			if i >= rows || j >= columns { return nil }
			return values[i][j]
		}
		set(newValue) {
			if i >= rows && j >= columns { return }
			if let newValue {
				values[i][j] = newValue
			}
		}
	}
	
	public enum MatrixError: Error {
		case nonMatchingDimensions
		case singularMatrix
	}
}
