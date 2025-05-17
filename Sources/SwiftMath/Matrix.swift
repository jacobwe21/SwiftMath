//
//  Matrix.swift
//  MySwift
//
//  Created by Jacob W Esselstyn on 1/19/25.
//
import Foundation

public struct Matrix<T: BinaryFloatingPoint>: Hashable, CustomStringConvertible {
	
	public let rows: Int
	public let columns: Int
	public private(set) var values: [[T]]
	
	/// Creates a matrix where the outer array is for rows and the inner array is for columns. Throws errors.
	public init(_ values: [[T]], mustBeSquare: Bool) throws {
		guard !values.isEmpty else {
			throw MatrixError.emptyMatrix
		}
		// Verify that the matrix is square
		if mustBeSquare {
			for row in values {
				guard row.count == values.count else {
					throw MatrixError.notSquareMatrix
				}
			}
		}
		self.rows = values.count
		self.columns = values.first!.count
		self.values = values
	}
	/// Creates a matrix where the outer array is for rows and the inner array is for columns.
	public init(_ values: [[T]]) {
		guard !values.isEmpty else {
			fatalError("Matrix Initialization Failure. Matrix is Empty.")
		}
		self.rows = values.count
		self.columns = values.first!.count
		self.values = values
	}
	/// Initialize a matrix of size n x m with a default value
	public init(rows n: Int, columns m: Int, defaultValue: T = 0) {
		self.values = Array(repeating: Array(repeating: defaultValue, count: m), count: n)
		self.rows = n
		self.columns = m
	}
	/// Initialize a matrix of size n x n with a default value
	public init(size n: Int, defaultValue: T = 0) {
		self.values = Array(repeating: Array(repeating: defaultValue, count: n), count: n)
		self.rows = n; self.columns = n
	}
	/// Return the identity matrix of the same size
	public static func identity(size n: Int) -> Matrix<T> {
		var result = Matrix<T>(size: n)
		for i in 0..<n {
			result[i, i] = 1 as! T
		}
		return result
	}
	
	/// Check if rows and columns of each matrix match.
	public static func areDimensionsEqual(lhs: Self, rhs: Self) -> Bool {
		lhs.rows == rhs.rows && lhs.columns == rhs.columns
	}
	
	/// Get a row as an array
	public func row(_ index: Int) -> [T] {
		guard index >= 0 && index < rows else {
			fatalError("Row index out of bounds: \(index)")
		}
		return values[index]
	}
	/// Get a column as an array
	public func column(_ index: Int) -> [T] {
		guard index >= 0 && index < columns else {
			fatalError("Column index out of bounds: \(index)")
		}
		return (0..<rows).map { values[$0][index] }
	}
	
	/// Elementwise addition
	public static func + (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = lhs.values[i][j] - rhs.values[i][j]
			}
		}
		return result
	}
	/// Elementwise subtraction
	public static func - (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard Matrix.areDimensionsEqual(lhs: lhs, rhs: rhs) else { throw MatrixError.nonMatchingDimensions }
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = lhs.values[i][j] - rhs.values[i][j]
			}
		}
		return result
	}
	
	/// Scalar multiplication
	public static func * (lhs: Matrix<T>, rhs: T) -> Matrix<T> {
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]*rhs
			}
		}
		return result
	}
	/// Scalar multiplication
	public static func * (lhs: T, rhs: Matrix<T>) -> Matrix<T> {
		var result = rhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]*lhs
			}
		}
		return result
	}
	
	/// Elementwise multiplication of two matricies
	public mutating func multiplyElementwise(by other: Self) {
		for i in 0..<rows {
			for j in 0..<columns {
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
	
	/// Perform Matrix Multiplication
	public static func * (lhs: Matrix<T>, rhs: Matrix<T>) throws -> Matrix<T> {
		guard lhs.columns == rhs.rows else { throw MatrixError.nonMatchingDimensions }
		var result = Matrix<T>(rows: lhs.rows, columns: rhs.columns)
		for i in 0..<lhs.rows {
			for j in 0..<rhs.columns {
				for k in 0..<lhs.columns {
					result.values[i][j] += lhs.values[i][k] * rhs.values[k][j]
				}
			}
		}
		return result
	}
	
	/// Returns the transpose of the matrix (swapping rows for columns)
	public func transpose() -> Self {
		guard rows > 0 && columns > 0 else { return self }
		var result = Matrix(rows: columns, columns: rows)
		for i in 0..<rows {
			for j in 0..<columns {
				result.values[j][i] = self.values[i][j]
			}
		}
		return result
	}
	
	/// Returns `true` if matrix is symmetric
	public var isSymmetric: Bool {
		if isSquare { return self == self.transpose() } else { return false }
	}
	
	/// Check if the matrix is square.
	private var isSquare: Bool { rows == columns }
	
	/// Determinant of the matrix (only defined for square matrices).
	public func determinant() throws -> T {
		guard isSquare else { throw MatrixError.notSquareMatrix }
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
		if n == 3 { // Base case for 3x3 matrix
			let a = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
			let b = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
			let c = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
			return a - b + c
		}
	
		var det: T = 0
		// Calculate determinant using row reduction (Gaussian elimination) - O(n^3)
		// TO-DO
		// Calculate determinant using cofactor expansion - O(n!)
		// Calculate determinant using LU decomposition?
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
	
	/// Check if indices are valid
	private func isValidIndex(row: Int, col: Int) -> Bool {
		return row >= 0 && row < rows && col >= 0 && col < columns
	}
	/// Fetches or updates a value in the matrix. Will crash if out-of-range.
	public subscript(i: Int, j: Int) -> T {
		get {
			guard isValidIndex(row: i, col: j) else {
				fatalError("Matrix index out of bounds: (\(i), \(j))")
			}
			return values[i][j]
		}
		set(newValue) {
			guard isValidIndex(row: i, col: j) else {
				fatalError("Matrix index out of bounds: (\(i), \(j))")
			}
			values[i][j] = newValue
		}
	}
//	/// Fetches or updates a value in the matrix. Returns `nil` if out-of-range, and does nothing when setting a value if out-of-range or if new value is `nil`.
//	public subscript(i: Int, j: Int) -> T? {
//		get {
//			guard isValidIndex(row: i, col: j) else { return nil }
//			return values[i][j]
//		}
//		set(newValue) {
//			guard isValidIndex(row: i, col: j) && newValue != nil else { return }
//			if let newValue {
//				values[i][j] = newValue
//			}
//		}
//	}
	
	// Error types for matrix operations
	public enum MatrixError: Error {
		case nonMatchingDimensions
		case singularMatrix
		case emptyMatrix
		case notSquareMatrix
		case indexOutOfRange
	}
	
	/// Convert matrix to a string representation
	public var description: String {
		var result = ""
		for i in 0..<rows {
			result += "["
			for j in 0..<columns {
				result += "\(values[i][j])"
				if j < columns - 1 {
					result += ", "
				}
			}
			result += "]\n"
		}
		return result
	}
}
