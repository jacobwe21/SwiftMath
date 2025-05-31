//
//  Matrix.swift
//  MySwift
//
//  Created by Jacob W Esselstyn on 1/19/25.
//
import Foundation
import simd
import Accelerate

public struct Matrix<T: BinaryFloatingPoint>: Equatable, CustomStringConvertible {
	
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
	/// Creates a matrix from a single array, either as a row or as a column.
	public init(_ values: [T], isRow: Bool) {
		guard !values.isEmpty else {
			fatalError("Matrix Initialization Failure. Matrix is Empty.")
		}
		if isRow {
			self.rows = 1
			self.columns = values.count
			self.values = [values]
		} else {
			self.rows = values.count
			self.columns = 1
			self.values = values.map({[$0]})
		}
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
	public init(_ simdMatrix: simd_double2x2) where T == Double {
		self.rows = 2
		self.columns = 2
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y],
			[m.columns.1.x, m.columns.1.y],
		]
	}
	public init(_ simdMatrix: simd_double3x3) where T == Double {
		self.rows = 3
		self.columns = 3
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y, m.columns.0.z],
			[m.columns.1.x, m.columns.1.y, m.columns.1.z],
			[m.columns.2.x, m.columns.2.y, m.columns.2.z],
		]
	}
	public init(_ simdMatrix: simd_double4x4) where T == Double {
		self.rows = 4
		self.columns = 4
		let m = simdMatrix.transpose
		self.values = [
			[m.columns.0.x, m.columns.0.y, m.columns.0.z, m.columns.0.w],
			[m.columns.1.x, m.columns.1.y, m.columns.1.z, m.columns.1.w],
			[m.columns.2.x, m.columns.2.y, m.columns.2.z, m.columns.2.w],
			[m.columns.3.x, m.columns.3.y, m.columns.3.z, m.columns.3.w],
		]
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
	
	/// Returns `true` if both matrices are approximately equivalent to 10^-4
	public static func ~= (lhs: Matrix<T>, rhs: Matrix<T>) -> Bool {
		guard lhs.sizeDescription == rhs.sizeDescription else { return false }
		let result = try! lhs - rhs
		for r in 0..<lhs.rows {
			for c in 0..<lhs.columns {
				if abs(result.values[r][c]) > 0.0001 {
					return false
				}
			}
		}
		return true
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
	/// Scalar division
	public static func / (lhs: Matrix<T>, rhs: T) -> Matrix<T> {
		var result = lhs
		for i in 0..<result.rows {
			for j in 0..<result.columns {
				result.values[i][j] = result.values[i][j]/rhs
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
	
	/// Returns the matrix as a single array, only if matrix is a single column or a single row.
	public func makeVectorArray() throws -> [T] {
		guard rows == 1 || columns == 1 else { throw MatrixError.nonMatchingDimensions }
		if rows == 1 {
			return values[0]
		} else {
			return self.transpose().values[0]
		}
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
	}
	
	public func cofactorMatrix() throws -> Matrix<T> {
		guard isSquare else {
			throw MatrixError.notSquareMatrix
		}
		let n = rows
		guard n > 1 else { throw MatrixError.indexOutOfRange }
		if n == 2 {
			return Matrix([
				[values[1][1], -values[1][0]],
				[-values[0][1], values[0][0]]
			])
		}
		var cofactorMatrix = [[T]](repeating: [T](repeating: 0, count: n), count: n)
		for i in 0..<n {
			for j in 0..<n {
				let minor = minorMatrix(removingRow: i, removingColumn: j)
				let sign: T = ((i + j) % 2 == 0) ? 1 : -1
				cofactorMatrix[i][j] = sign * Matrix.calculateDeterminant(for: minor)
			}
		}
		return Matrix(cofactorMatrix)
	}
	public func adjoint() throws -> Matrix<T> {
		// Transpose the cofactor matrix to get the adjoint
		return try self.cofactorMatrix().transpose()
	}
	private func minorMatrix(removingRow rowToRemove: Int, removingColumn columnToRemove: Int) -> [[T]] {
		values.enumerated().compactMap { i, row in
			guard i != rowToRemove else { return nil }
			let filteredRow = row.enumerated().compactMap { j, value in
				j != columnToRemove ? value : nil
			}
			return filteredRow
		}
	}
	
	public func removing(rows rowsToRemove: [Int], columns columnsToRemove: [Int]) throws -> Matrix<T> {
		for row in rowsToRemove {
			if row < 0 || row >= rows { throw MatrixError.indexOutOfRange }
		}
		for column in columnsToRemove {
			if column < 0 || column >= columns { throw MatrixError.indexOutOfRange }
		}
		let rowSet = Set(rowsToRemove)
		let columnSet = Set(columnsToRemove)
		let filteredValues: [[T]] = values.enumerated().compactMap { (i, row) in
			guard !rowSet.contains(i) else { return nil }
			let filteredRow = row.enumerated().compactMap { (j, value) in
				columnSet.contains(j) ? nil : value
			}
			return filteredRow
		}
		return Matrix(filteredValues)
	}
	
	/// Determinant of the matrix (only defined for square matrices).
	public func inverse() throws -> Matrix<T> {
		guard isSquare else { throw MatrixError.notSquareMatrix }
		
		let det = Matrix.calculateDeterminant(for: values)
		if det == 0 { throw MatrixError.singularMatrix }
		
		let size = rows
		if size == 1 { return Matrix([[(1/values[0][0])]]) }
		if size == 3 {
			let m = values

			let a = m[0][0], b = m[0][1], c = m[0][2]
			let d = m[1][0], e = m[1][1], f = m[1][2]
			let g = m[2][0], h = m[2][1], i = m[2][2]

			let A = e * i - f * h
			let B = -(d * i - f * g)
			let C = d * h - e * g
			let D = -(b * i - c * h)
			let E = a * i - c * g
			let F = -(a * h - b * g)
			let G = b * f - c * e
			let H = -(a * f - c * d)
			let I = a * e - b * d
			
			let det = a * A + b * B + c * C
			guard det != 0 else {
				throw MatrixError.singularMatrix
			}
			let resultValues: [[T]] = [
				[A, D, G],
				[B, E, H],
				[C, F, I]
			]
			return Matrix(resultValues)/det
		}
		if size < 5 {
			return (1/det)*(try! self.adjoint())
		}
		
		var augmented = values

		// Append identity matrix to the right of the original
		for i in 0..<size {
			augmented[i] += (0..<size).map { $0 == i ? 1 : 0 }.map(T.init)
		}

		// Perform Gauss-Jordan elimination
		for i in 0..<size {
			// Find pivot
			var pivotRow = i
			for j in i+1..<size where abs(augmented[j][i]) > abs(augmented[pivotRow][i]) {
				pivotRow = j
			}

			// If pivot is zero, matrix is singular
			if augmented[pivotRow][i] == 0 {
				throw MatrixError.singularMatrix
			}

			// Swap rows if needed
			if i != pivotRow {
				augmented.swapAt(i, pivotRow)
			}

			// Normalize pivot row
			let pivot = augmented[i][i]
			for j in 0..<2 * size {
				augmented[i][j] /= pivot
			}

			// Eliminate other rows
			for k in 0..<size where k != i {
				let factor = augmented[k][i]
				for j in 0..<2 * size {
					augmented[k][j] -= factor * augmented[i][j]
				}
			}
		}

		// Extract right half as the inverse
		let inverseValues = augmented.map { Array($0[size..<(2*size)]) }
		return Matrix(inverseValues)
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
		// TO-DO
		// Calculate determinant using row reduction (Gaussian elimination) - O(n^3)
		//return rowEcholonForm?.mainDiagonal?.reduce(1, *) ?? 0
		// Calculate determinant using LU decomposition?
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
	public enum MatrixError: String, Error {
		case nonMatchingDimensions = "Dimensions of matrices do not match."
		case singularMatrix = "Matrix is singular."
		case emptyMatrix = "Matrix is empty."
		case notSquareMatrix = "Matrix is not square, and must be square for this operation."
		case indexOutOfRange = "Index is out-of-range."
	}
	
	/// Display row count x column count
	public var sizeDescription: String {
		"\(rows)x\(columns)"
	}
	
	/// Convert matrix to a string representation
	public var description: String {
		var result = "\(rows)x\(columns)\n"
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

// MARK: Large Sparse Matrices
public extension SparseMatrix_Double {
	
	enum MatrixUpdateError: Error {
		case couldNotFindData
		case invalidRowIndex
		case invalidColumnIndex
	}
	mutating func updateMatrix(at row: Int, column: Int, with newValue: Double) throws {
		guard row >= 0 else { throw MatrixUpdateError.invalidRowIndex }
		guard column >= 0 else { throw MatrixUpdateError.invalidColumnIndex }
		guard row < Int(self.structure.rowCount) else { throw MatrixUpdateError.invalidRowIndex }
		guard column < Int(self.structure.columnCount) else { throw MatrixUpdateError.invalidColumnIndex }
		
		// loop through the matrix to find the index for the data to update
		let nnz = self.structure.columnStarts[Int(self.structure.columnCount)]
		var i = 0
		var currentColumn: Int = 0
		var nextColStarts = self.structure.columnStarts[1]
		while i < nnz {
			if i == nextColStarts {
				currentColumn += 1
				if currentColumn > column { throw MatrixUpdateError.couldNotFindData }
				nextColStarts = self.structure.columnStarts[currentColumn + 1]
			}
			let rowIndex = Int(self.structure.rowIndices[i])
			if rowIndex == row && currentColumn == column {
				self.data[i] = newValue
				return
			}
			i += 1
		}
		throw MatrixUpdateError.couldNotFindData
	}
	
	// https://stackoverflow.com/questions/53006935/how-to-read-values-of-sparsematrix-double-in-swift-4
	var array: [[Double]] {

		let rows = Int(self.structure.rowCount)
		let columns = Int(self.structure.columnCount)
		let nnz = self.structure.columnStarts[Int(self.structure.columnCount)]
				
		var M = Array(repeating: Array(repeating: 0.0, count: columns), count: rows)

		var i = 0
		var currentColumn: Int = 0
		var nextColStarts = self.structure.columnStarts[1]
		while i < nnz {
			if i == nextColStarts {
				currentColumn += 1
				nextColStarts = self.structure.columnStarts[currentColumn + 1]
			}

			let rowIndex = Int(self.structure.rowIndices[i])
			M[rowIndex][currentColumn] = self.data[i]
			i += 1
		}
		return M
	}
	
	subscript (row: Int, column: Int) -> Double {
		get {
			self.array[row][column]
		} set {
			do {
				try self.updateMatrix(at: row, column: column, with: newValue)
			} catch {
				fatalError(error.localizedDescription)
			}
		}
	}
}
//public struct SparseMatrix<T: FloatingPoint> {
//	public var columnCount: Int32
//	public var rowCount: Int32
//	public var columnStarts: [Int]
//	public var rowIndices: [Int32]
//	public var values: [T]
//	public var structure: SparseMatrixStructure
//
//	public init(_ values: [[T]]) {
//
//	}
//
//	public init(rowCount: Int32, columnCount: Int32, columnStarts: [Int], rowIndices: [Int32], values: [T], symmetric: Bool, upperTriangle: Bool = true) {
//		self.columnCount = columnCount
//		self.rowCount = rowCount
//		self.columnStarts = columnStarts
//		self.rowIndices = rowIndices
//		self.values = values
//
//		var attributes = SparseAttributes_t()
//		if symmetric {
//			if upperTriangle {
//				attributes.triangle = SparseUpperTriangle
//			} else {
//				attributes.triangle = SparseLowerTriangle
//			}
//			attributes.kind = SparseSymmetric
//		}
//		self.columnStarts = columnStarts
//		self.rowIndices = rowIndices
//		structure = self.rowIndices.withUnsafeMutableBufferPointer { rowIndicesPtr in
//			self.columnStarts.withUnsafeMutableBufferPointer { columnStartsPtr in
//				return SparseMatrixStructure(
//					rowCount: rowCount,
//					columnCount: columnCount,
//					columnStarts: columnStartsPtr.baseAddress!,
//					rowIndices: rowIndicesPtr.baseAddress!,
//					attributes: attributes,
//					blockSize: 1
//				)
//			}
//		}
//	}
//}

// MARK: System Solving - Test Code Fails
//	mutating func runAnalysis() throws {
//		results = nil
//		let globalMatrixSize: Int32 = nodeCount*3
////		// Organize entries for the global stiffness matrix. Evaluate element stiffness coefficients for each element, and assign subcript/index (i = Force ID (row), j = Degree of Freedom ID (column))
////		var memberStiffnessEntries: [SIMD2<Int32>:Double] = [:] // Key is (i,j), Value is Stiffness
////		for member in members {
//////			var memberStiffnessValues = member.values
//////			memberStiffnessValues.withUnsafeMutableBufferPointer { valuesPtr in
//////				let M = SparseMatrix_Double(structure: member.structureOfMatrix, data: valuesPtr.baseAddress!)
//////				SparseCleanup(M)
//////				//memberStiffnessEntries.append()
//////			}
////			let memberStiffnessMatrix = member.stiffnessMatrix
////			for i in memberStiffnessMatrix.indices {
////				for j in memberStiffnessMatrix[i].indices {
////					let iNode = i<3 ? member.node1 : member.node2
////					let iNodeIndex = nodes.firstIndex(of: iNode)!
////					let forceID: Int32 = Int32(iNodeIndex*3+i%3)
////					let jNode = j<3 ? member.node1 : member.node2
////					let jNodeIndex = nodes.firstIndex(of: jNode)!
////					let dofID: Int32 = Int32(jNodeIndex*3+j%3)
////					let key = SIMD2(x: forceID, y: dofID)
////					if let currentValue = memberStiffnessEntries[key] {
////						memberStiffnessEntries.updateValue(currentValue+memberStiffnessMatrix[i][j], forKey: key)
////					} else {
////						memberStiffnessEntries.updateValue(memberStiffnessMatrix[i][j], forKey: key)
////					}
////				}
////			}
////		}
////
////		// Create the global stiffness matrix based on the member stiffnesses.
////		var row: [Int32] = []
////		var column: [Int32] = []
////		var values: [Double] = []
////		for entry in memberStiffnessEntries {
////			row.append(entry.key.x)
////			column.append(entry.key.y)
////			values.append(entry.value)
////		}
////		var attributes = SparseAttributes_t()
////		let blockCount = values.count
////		let blockSize = 1
////		// `SparseConvertFromCoordinate` sums all duplicate entries, assembling the matrix correctly.
////		let K = SparseConvertFromCoordinate(globalMatrixSize, globalMatrixSize,
////											blockCount, UInt8(blockSize),
////											attributes,
////											&row, &column,
////											&values)
//
//
//		// Create Kff matrix for displacements
//		var kffMatrixSize: Int32 = globalMatrixSize
//		for node in nodes {
//			if node.supportResistance.forFx { kffMatrixSize -= 1 }
//			if node.supportResistance.forFy { kffMatrixSize -= 1 }
//			if node.supportResistance.forMz { kffMatrixSize -= 1 }
//		}
//		if kffMatrixSize == globalMatrixSize {
//			print("Under Constrained. No supports. Singluar Matrix.")
//			throw AnalysisError.singularMatrix
//		}
//		if kffMatrixSize == 0 {
//			print("Over constrained, stiffness matrix is empty")
//			throw AnalysisError.overconstrained
//		}
//		var stiffnessesForKff: [SIMD2<Int32>:Double] = [:]
//		var forcesVectorF: [Int32:Double] = [:]
//		var kffNodeDisplacementDictionary: [Int32:SIMD2<Double>] = [:]
//		var iSkip: Int = 0
//		var jSkip: Int = 0
//		for member in members {
//			let memberStiffnessMatrix = member.stiffnessMatrix()
//			for i in 0..<memberStiffnessMatrix.rows {
//				if i == 0 { if member.node1.supportResistance.forFx { iSkip += 1; continue }}
//				else if i == 1 { if member.node1.supportResistance.forFy { iSkip += 1; continue }}
//				else if i == 2 { if member.node1.supportResistance.forMz { iSkip += 1; continue }}
//				else if i == 3 { if member.node2.supportResistance.forFx { iSkip += 1; continue }}
//				else if i == 4 { if member.node2.supportResistance.forFy { iSkip += 1; continue }}
//				else if i == 5 { if member.node2.supportResistance.forMz { iSkip += 1; continue }}
//				else { throw AnalysisError.indexOutOfRange }
//				for j in 0..<memberStiffnessMatrix.columns {
//					if j == 0 { if member.node1.supportResistance.forFx { jSkip += 1; continue }}
//					else if j == 1 { if member.node1.supportResistance.forFy { jSkip += 1; continue }}
//					else if j == 2 { if member.node1.supportResistance.forMz { jSkip += 1; continue }}
//					else if j == 3 { if member.node2.supportResistance.forFx { jSkip += 1; continue }}
//					else if j == 4 { if member.node2.supportResistance.forFy { jSkip += 1; continue }}
//					else if j == 5 { if member.node2.supportResistance.forMz { jSkip += 1; continue }}
//					else { throw AnalysisError.indexOutOfRange }
//					let iNode = i<3 ? member.node1 : member.node2
//					let iNodeIndex = nodes.firstIndex(of: iNode)!
//					let forceID: Int32 = Int32(iNodeIndex*3+i%3-iSkip)
//					let jNode = j<3 ? member.node1 : member.node2
//					let jNodeIndex = nodes.firstIndex(of: jNode)!
//					let dofID: Int32 = Int32(jNodeIndex*3+j%3-jSkip)
//
//					let key = SIMD2(x: forceID, y: dofID)
//					if let currentValue = stiffnessesForKff[key] {
//						stiffnessesForKff.updateValue(currentValue+memberStiffnessMatrix[i,j], forKey: key)
//					} else {
//						stiffnessesForKff.updateValue(memberStiffnessMatrix[i, j], forKey: key)
//					}
//
//					let forceAtNode: Double
//					if i == 0 { forceAtNode = member.loads1.x }
//					else if i == 1 { forceAtNode = member.loads1.y }
//					else if i == 2 { forceAtNode = member.loads1.z }
//					else if i == 3 { forceAtNode = member.loads2.x }
//					else if i == 4 { forceAtNode = member.loads2.y }
//					else if i == 5 { forceAtNode = member.loads2.z }
//					else { throw AnalysisError.indexOutOfRange }
//					if let currentValue = forcesVectorF[forceID] {
//						forcesVectorF.updateValue(currentValue+forceAtNode, forKey: forceID)
//					} else {
//						forcesVectorF.updateValue(forceAtNode, forKey: forceID)
//					}
//
//					if j%3 == 0 {
//						kffNodeDisplacementDictionary.updateValue(jNode.id, forKey: dofID)
//					}
//				}
//			}
//		}
//		var rowKff: [Int32] = []
//		var columnKff: [Int32] = []
//		var valuesKff: [Double] = []
//		for entry in stiffnessesForKff {
//			rowKff.append(entry.key.x)
//			columnKff.append(entry.key.y)
//			valuesKff.append(entry.value)
//		}
//		let attributesKff = SparseAttributes_t()
//		let blockCountKff = valuesKff.count
//		let blockSizeKff = 1
//		// `SparseConvertFromCoordinate` sums all duplicate entries, assembling the matrix correctly.
//		let Kff = SparseConvertFromCoordinate(kffMatrixSize, kffMatrixSize,
//											blockCountKff, UInt8(blockSizeKff),
//											attributesKff,
//											&rowKff, &columnKff,
//											&valuesKff)
//		let factoredKff = SparseFactor(SparseFactorizationCholesky, Kff)
//		SparseCleanup(Kff)
//		SparseCleanup(factoredKff)
//
//		// Create Ksf matrix for reaction forces
//		var stiffnessesForKsf: [SIMD2<Int32>:Double] = [:]
//		var reactionOffsetVector: [Int32:Double] = [:]
//		var ksfNodeReactionDictionary: [Int32:SIMD2<Double>] = [:]
//		iSkip = 0
//		jSkip = 0
//		for member in members {
//			let memberStiffnessMatrix = member.stiffnessMatrix()
//			for i in 0..<memberStiffnessMatrix.rows {
//				if i == 0 { if !member.node1.supportResistance.forFx { iSkip += 1; continue }}
//				else if i == 1 { if !member.node1.supportResistance.forFy { iSkip += 1; continue }}
//				else if i == 2 { if !member.node1.supportResistance.forMz { iSkip += 1; continue }}
//				else if i == 3 { if !member.node2.supportResistance.forFx { iSkip += 1; continue }}
//				else if i == 4 { if !member.node2.supportResistance.forFy { iSkip += 1; continue }}
//				else if i == 5 { if !member.node2.supportResistance.forMz { iSkip += 1; continue }}
//				else { throw AnalysisError.indexOutOfRange }
//				for j in 0..<memberStiffnessMatrix.columns {
//					if j == 0 { if !member.node1.supportResistance.forFx { jSkip += 1; continue }}
//					else if j == 1 { if !member.node1.supportResistance.forFy { jSkip += 1; continue }}
//					else if j == 2 { if !member.node1.supportResistance.forMz { jSkip += 1; continue }}
//					else if j == 3 { if !member.node2.supportResistance.forFx { jSkip += 1; continue }}
//					else if j == 4 { if !member.node2.supportResistance.forFy { jSkip += 1; continue }}
//					else if j == 5 { if !member.node2.supportResistance.forMz { jSkip += 1; continue }}
//					else { throw AnalysisError.indexOutOfRange }
//					let iNode = i<3 ? member.node1 : member.node2
//					let iNodeIndex = nodes.firstIndex(of: iNode)!
//					let forceID: Int32 = Int32(iNodeIndex*3+i%3-iSkip)
//					let jNode = j<3 ? member.node1 : member.node2
//					let jNodeIndex = nodes.firstIndex(of: jNode)!
//					let dofID: Int32 = Int32(jNodeIndex*3+j%3-jSkip)
//
//					let key = SIMD2(x: forceID, y: dofID)
//					if let currentValue = stiffnessesForKsf[key] {
//						stiffnessesForKsf.updateValue(currentValue+memberStiffnessMatrix[i, j], forKey: key)
//					} else {
//						stiffnessesForKsf.updateValue(memberStiffnessMatrix[i, j], forKey: key)
//					}
//
//					let forceAtNode: Double
//					if i == 0 { forceAtNode = member.loads1.x }
//					else if i == 1 { forceAtNode = member.loads1.y }
//					else if i == 2 { forceAtNode = member.loads1.z }
//					else if i == 3 { forceAtNode = member.loads2.x }
//					else if i == 4 { forceAtNode = member.loads2.y }
//					else if i == 5 { forceAtNode = member.loads2.z }
//					else { throw AnalysisError.indexOutOfRange }
//					if let currentValue = forcesVectorF[forceID] {
//						reactionOffsetVector.updateValue(currentValue+forceAtNode, forKey: forceID)
//					} else {
//						reactionOffsetVector.updateValue(forceAtNode, forKey: forceID)
//					}
//
//					if i%3 == 0 {
//						ksfNodeReactionDictionary.updateValue(iNode.id, forKey: forceID)
//					}
//				}
//			}
//		}
//		// Solve for reactions using Ksf matrix and displacements vector
//		let ksfMatrixSize = globalMatrixSize - kffMatrixSize
//		var rowKsf: [Int32] = []
//		var columnKsf: [Int32] = []
//		var valuesKsf: [Double] = []
//		for entry in stiffnessesForKsf {
//			rowKsf.append(entry.key.x)
//			columnKsf.append(entry.key.y)
//			valuesKsf.append(entry.value)
//		}
//		let attributesKsf = SparseAttributes_t()
//		let blockCountKsf = valuesKsf.count
//		let blockSizeKsf = 1
//		let Ksf = SparseConvertFromCoordinate(ksfMatrixSize, ksfMatrixSize,
//											blockCountKsf, UInt8(blockSizeKsf),
//											attributesKsf,
//											&rowKsf, &columnKsf,
//											&valuesKsf)
//		SparseCleanup(Ksf)
//
//		// Perform matrix operations using Kff and Ksf
//		guard forcesVectorF.count == kffMatrixSize else {
//			print("forcesVectorF != kffMatrixSize"); throw AnalysisError.matrixSizeMisMatch
//		}
//		guard reactionOffsetVector.count == ksfMatrixSize else {
//			print("forcesVectorF != kffMatrixSize"); throw AnalysisError.matrixSizeMisMatch
//		}
//		var forcesVectorFValues: [Double] = Array(forcesVectorF.values)
//		let reactionOffsetVectorValues: [Double] = Array(reactionOffsetVector.values)
//		var displacements: [Double] = Array(repeating: 0, count: Int(kffMatrixSize))
//		var reactions: [Double] = Array(repeating: 0, count: reactionOffsetVectorValues.count)
//		forcesVectorFValues.withUnsafeMutableBufferPointer { bPtr in
//			// Forces
//			let b = DenseVector_Double(count: kffMatrixSize,
//											data: bPtr.baseAddress!)
//			displacements.withUnsafeMutableBufferPointer { xPtr in
//				// Displacements
//				let x = DenseVector_Double(count: kffMatrixSize, data: xPtr.baseAddress!)
//				SparseSolve(factoredKff, b, x) // Solve for displacements
//
//				reactions.withUnsafeMutableBufferPointer { xPtr in
//					let y: DenseVector_Double = DenseVector_Double(count: ksfMatrixSize, data: xPtr.baseAddress!)
//					SparseMultiply(Ksf,x,y) // Multiply for reactions
//				}
//			}
//		}
//		// Adjust reactions
//		for i in reactions.indices {
//			reactions[i] = reactions[i] + reactionOffsetVectorValues[i]
//		}
//		SparseCleanup(Kff)
//		SparseCleanup(Ksf)
//
//		var nodeDisplacements: [SIMD2<Double>:SIMD3<Double>] = [:]
//		for i in displacements.indices {
//			if let item = kffNodeDisplacementDictionary[Int32(i)] {
//				let d = SIMD3(displacements[i], displacements[i+1], displacements[i+2])
//				nodeDisplacements.updateValue(d, forKey: item)
//			}
//		}
//		var nodeReactions: [SIMD2<Double>:SIMD3<Double>] = [:]
//		for i in reactions.indices {
//			if let item = ksfNodeReactionDictionary[Int32(i)] {
//				let r = SIMD3(reactions[i], reactions[i+1], reactions[i+2])
//				nodeReactions.updateValue(r, forKey: item)
//			}
//		}
//		results = Results(displacements: nodeDisplacements, reactions: nodeReactions)
//	}
