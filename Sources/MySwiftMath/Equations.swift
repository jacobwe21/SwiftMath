//
//  Equations.swift
//
//  Created by Jacob W Esselstyn on 12/29/22.
//

import SwiftUI
import Accelerate

public protocol MathEquation: CustomStringConvertible, Sendable {
	func callAsFunction(_ x: Double) -> Double
	func makeDerivative() -> MathEquation
	func integrate(plus c: Double) -> MathEquation
	mutating func scale(by rhs: Double)
	func scaled(by rhs: Double) -> Self
	func criticalPoints(in range: ClosedRange<Double>?) -> [Double]
	func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double) -> (min: Double, max: Double)
}
extension MathEquation {
	public func definiteIntegral(over domain: ClosedRange<Double>, absoluteTolerance: Double = 1.0e-8, relativeTolerance: Double = 1.0e-2) throws -> Double {
		let q = Quadrature(integrator: .qags(maxIntervals: 12), absoluteTolerance: absoluteTolerance, relativeTolerance: relativeTolerance)
		let result = q.integrate(over: domain, integrand: callAsFunction)
		switch result {
		case .success((let integralResult, let estimatedAbsoluteError)):
			print("quadrature success:", integralResult,
				  estimatedAbsoluteError)
			return integralResult
		case .failure(let error):
			print("quadrature error:", error.errorDescription)
			throw error
		}
	}
}

public struct Math {
	
	public struct MultiEQ: MathEquation {
		
		public init(segments: [Segment]) {
			self.segments = segments
		}
		public init(equations: [any MathEquation]) {
			var newSegments = [Segment]()
			for eq in equations {
				newSegments.append(Segment(eq: eq))
			}
			self.segments = newSegments
		}
		
		public struct Segment: Sendable {
			var eq: any MathEquation
			let xStart: Double
			let xEnd: Double
			let xStartIsInclusive: Bool
			let xEndIsInclusive: Bool
			init(eq: any MathEquation, xStart: Double = -Double.infinity, xEnd: Double = Double.infinity, xStartIsInclusive: Bool = false, xEndIsInclusive: Bool = false) {
				self.eq = eq
				self.xStart = xStart
				self.xEnd = xEnd
				self.xStartIsInclusive = xStartIsInclusive
				self.xEndIsInclusive = xEndIsInclusive
			}
			public init(equation: any MathEquation, xStart: Double = -Double.infinity, xEnd: Double = Double.infinity, xStartIsInclusive: Bool, xEndIsInclusive: Bool) {
				self.init(eq: equation, xStart: xStart, xEnd: xEnd, xStartIsInclusive: xStartIsInclusive, xEndIsInclusive: xEndIsInclusive)
			}
			public init(equation: any MathEquation, range: ClosedRange<Double>) {
				self.init(eq: equation, xStart: range.lowerBound, xEnd: range.upperBound, xStartIsInclusive: true, xEndIsInclusive: true)
			}
			public mutating func scale(by rhs: Double) {
				eq.scale(by: rhs)
			}
		}
		public private(set) var segments: [Segment]
		
		public func callAsFunction(_ x: Double) -> Double {
			var result = 0.0
			for s in segments {
				if s.xStart == x && !s.xStartIsInclusive { continue }
				if s.xEnd == x && !s.xEndIsInclusive { continue }
				if s.xStart <= x && s.xEnd >= x { result += s.eq(x) }
			}
			return result
		}
		
		public func makeDerivative() -> any MathEquation {
			var newSegments: [Segment] = []
			for s in segments {
				if s.eq is Math.Impulse { continue }
				if s.eq is Math.NullEQ { continue }
				var segment = s
				segment.eq = s.eq.makeDerivative()
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		
		/// Integrates simple equations faster and results in less additional segments.
		func quickIntegration(plus c: Double) -> (any MathEquation)? {
			if segments.count == 0 { return NullEQ().integrate(plus: c) }
			if segments.count == 1 {
				var tempSegment = segments.first!
				if tempSegment.eq is Math.Impulse {
					tempSegment = Segment(eq: Math.PolynomialEQ(y: tempSegment.eq(0)), xStart: tempSegment.xStart, xEnd: Double.infinity, xStartIsInclusive: true, xEndIsInclusive: false)
				} else {
					let newEQ = tempSegment.eq.integrate(plus: c)
					tempSegment.eq = newEQ
				}
				return Math.MultiEQ(segments: [tempSegment])
			}
			return nil
		}
		func intergrateImpulseSegment(s: Segment) -> Segment {
			guard s.eq is Math.Impulse else { fatalError("Segment \(s.eq) is not an Impulse") }
			if s.eq(0).isInfinite {
				return Segment(eq: Math.PolynomialEQ(y: 1), xStart: s.xStart, xEnd: Double.infinity, xStartIsInclusive: true, xEndIsInclusive: false)
			} else {
				return Segment(eq: Math.PolynomialEQ(y: s.eq(0)), xStart: s.xStart, xEnd: Double.infinity, xStartIsInclusive: true, xEndIsInclusive: false)
			}
		}
		
		public func integrateAtZero(plus c: Double, xStart: Double = -Double.infinity) -> any MathEquation {
			if let result = quickIntegration(plus: c) { return result }
			var newSegments: [Segment] = []
			let sortedSegments = segments.sorted(by: {
				if $0.xStart == $1.xStart {
					return $0.xStartIsInclusive
				} else {
					return $0.xStart<$1.xStart
				}
			})
			// Intergrate segments
			for s in sortedSegments {
				if s.eq is Math.Impulse {
					newSegments.append(intergrateImpulseSegment(s: s))
				} else {
					var segment = s
					segment.eq = s.eq.integrate(plus: 0)
					segment.eq = s.eq.integrate(plus: -segment.eq(0))
					newSegments.append(segment)
				}
			}
			// Adjust by c
			if c != 0 {
				let segment = Segment(eq: Math.PolynomialEQ(y: c), xStart: xStart, xEnd: Double.infinity, xStartIsInclusive: false, xEndIsInclusive: false)
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		public func integrateSmoothly(plus c: Double, xStart: Double = -Double.infinity) -> any MathEquation {
			if let result = quickIntegration(plus: c) { return result }
			var newSegments: [Segment] = []
			let sortedSegments = segments.sorted(by: {
				if $0.xStart == $1.xStart {
					return $0.xStartIsInclusive
				} else {
					return $0.xStart<$1.xStart
				}
			})
			// Intergrate segments
			for s in sortedSegments {
				if s.eq is Math.Impulse {
					newSegments.append(intergrateImpulseSegment(s: s))
				} else {
					var segment = s
					segment.eq = s.eq.integrate(plus: 0)
					segment.eq = s.eq.integrate(plus: -segment.eq(s.xStart))
					newSegments.append(segment)
					if s.xEnd.isFinite || s.xEnd.isInfinite && s.xEnd < 0 {
						let continuitySegment = Segment(eq: Math.PolynomialEQ(y: segment.eq(s.xEnd)), xStart: s.xEnd, xEnd: Double.infinity, xStartIsInclusive: !segment.xEndIsInclusive, xEndIsInclusive: false)
						newSegments.append(continuitySegment)
					}
				}
			}
			// Adjust by c
			if c != 0 {
				let segment = Segment(eq: Math.PolynomialEQ(y: c), xStart: xStart, xEnd: Double.infinity, xStartIsInclusive: false, xEndIsInclusive: false)
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		
		public func integrate(plus c: Double) -> any MathEquation {
			if let result = quickIntegration(plus: c) { return result }
			var newSegments: [Segment] = []
			let sortedSegments = segments.sorted(by: {
				if $0.xStart == $1.xStart {
					return $0.xStartIsInclusive
				} else {
					return $0.xStart<$1.xStart
				}
			})
			// Intergrate segments
			for s in sortedSegments {
				if s.eq is Math.Impulse {
					if s.eq(0).isInfinite {
						newSegments.append(intergrateImpulseSegment(s: s))
					} else {
						let segment = Segment(eq: Math.PolynomialEQ(y: s.eq(0)), xStart: s.xStart, xEnd: Double.infinity, xStartIsInclusive: true, xEndIsInclusive: false)
						newSegments.append(segment)
					}
				} else {
					var segment = s
					segment.eq = s.eq.integrate(plus: 0)
					let cForSegment = -(segment.eq(s.xStart)-segment.eq(0))
					segment.eq = s.eq.integrate(plus: cForSegment)
					newSegments.append(segment)
				}
			}
			// Adjust by c
			if c != 0 {
				let segment = Segment(eq: Math.PolynomialEQ(y: c), xStart: -Double.infinity, xEnd: Double.infinity, xStartIsInclusive: false, xEndIsInclusive: false)
				newSegments.append(segment)
			}
			return MultiEQ(segments: newSegments)
		}
		public mutating func scale(by rhs: Double) {
			for i in segments.indices {
				segments[i].eq.scale(by: rhs)
			}
		}
		public func scaled(by rhs: Double) -> Self {
			var newSegments: [Segment] = []
			for segment in segments {
				var newSegment = segment
				newSegment.eq.scale(by: rhs)
				newSegments.append(newSegment)
			}
			return Self.init(segments: newSegments)
		}
		
		public var minX: Double {
			return segments.min(by: {$0.xStart < $1.xStart})!.xStart
		}
		public var maxX: Double {
			return segments.max(by: {$0.xEnd < $1.xEnd})!.xEnd
		}
		public var description: String {
			var str = "{\n"
			for s in segments {
				if s.xStartIsInclusive && s.xEndIsInclusive {
					str += s.eq.description + " [\(s.xStart),\(s.xEnd)] \n"
				} else if !s.xStartIsInclusive && s.xEndIsInclusive {
					str += s.eq.description + " (\(s.xStart),\(s.xEnd)] \n"
				} else if s.xStartIsInclusive && !s.xEndIsInclusive {
					str += s.eq.description + " [\(s.xStart),\(s.xEnd)) \n"
				} else {
					str += s.eq.description + " (\(s.xStart),\(s.xEnd)) \n"
				}
			}
			return str+"}"
		}
		
		public func criticalPoints(in range: ClosedRange<Double>?) -> [Double] {
			var points = Set<Double>()
			for s in segments {
				if s.xStart.isFinite {
					points.insert(s.xStart)
					points.insert(s.xStart+Double.leastNonzeroMagnitude)
					points.insert(s.xStart-Double.leastNonzeroMagnitude)
				}
				if s.xEnd.isFinite {
					points.insert(s.xEnd)
					points.insert(s.xEnd+Double.leastNonzeroMagnitude)
					points.insert(s.xEnd-Double.leastNonzeroMagnitude)
				}
				points.formUnion(s.eq.criticalPoints(in: range))
			}
			if let range {
				points = points.filter { return range.contains($0) }
			}
			return points.map(\.self)
		}
		public func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double = 120) -> (min: Double, max: Double) {
			guard segments.count > 0 else { return (0, 0) }
			if segments.count == 1 {
				return segments.first!.eq.minmaxValues(in: range, num_dx: num_dx)
			}
			let xStart = max(range?.lowerBound ?? minX, minX)
			let xEnd = min(range?.upperBound ?? maxX, maxX)
			var minValue = callAsFunction(xStart)
			var maxValue = minValue
			for i in 0...Int(num_dx) {
				let z = Double(i)/num_dx // Percent
				if z > 1 || z < 0 { fatalError("z out of range") }
				let x = xStart + (xEnd - xStart)*z
				minValue = min(minValue, callAsFunction(x))
				maxValue = max(maxValue, callAsFunction(x))
			}
			for point in criticalPoints(in: range) {
				minValue = min(minValue, callAsFunction(point))
				maxValue = max(maxValue, callAsFunction(point))
			}
			return (minValue, maxValue)
		}
		
		/// Returns the Dirac Delta function.
		public static func diracDeltaEQ(x: Double) -> MultiEQ {
			MultiEQ(segments: [.init(eq: Impulse(term: Double.infinity), xStart: x, xEnd: x, xStartIsInclusive: true, xEndIsInclusive: true)])
		}
	}
	
	public struct PolynomialEQ: MathEquation, AdditiveArithmetic {
		public private(set) var terms: [Term]
		
		/// A straight horizontal line
		public init(y: Double) {
			self.terms = [Term(y, exp: 0)]
		}
		public init(terms: [Term]) {
			self.terms = terms
			organizeTerms()
		}
		public init(terms: Term...) {
			self.terms = terms
			organizeTerms()
		}
//		public init(_ str1: String) {
//			var str2 = str1.trimmingCharacters(in: CharacterSet(arrayLiteral: " "))
//			str2 = str2.replacingOccurrences(of: "-", with: "+-")
//			var termStrs = str2.components(separatedBy: CharacterSet(arrayLiteral: "+"))
//			termStrs.updateEach { s in
//				s.trimmingCharacters(in: CharacterSet(charactersIn: " "))
//			}
//			termStrs.removeAll { $0.isEmpty }
//			var terms = [Term]()
//			for t in termStrs {
//				let t2 = t.replacingOccurrences(of: "x^", with: "x")
//				let components = t2.components(separatedBy: CharacterSet(arrayLiteral: "x"))
//				if components.count == 2 {
//					let a = components.first!
//					let n = components.last!
//					terms.append(Term(Double(a) ?? 1, xToTheDouble: Double(n) ?? 1))
//				} else if components.count == 1 {
//					let a = components.first!
//					terms.append(Term(Double(a) ?? 1, xToThe: 0))
//				}
//			}
//			self.terms = terms
//		}
//		public init(_ polynomial: String) {
//			let pattern = #"([-+]?[0-9]*\.?[0-9]*)x\^?([0-9]*)"#
//			let regex = try! NSRegularExpression(pattern: pattern)
//			let matches = regex.matches(in: polynomial, range: NSRange(polynomial.startIndex..., in: polynomial))
//
//			var terms = [Term]()
//			for match in matches {
//				let coefficient = Double((Range(match.range(at: 1), in: polynomial).flatMap { String(polynomial[$0]) }) ?? "1") ?? 1.0
//				let degree = Double((Range(match.range(at: 2), in: polynomial).flatMap { String(polynomial[$0]) }) ?? "1") ?? 1.0
//				terms.append(Term(coefficient, xToTheDouble: degree))
//			}
//			self.terms = terms
//		}
		
		public func callAsFunction(_ x: Double) -> Double {
			var result = 0.0
			for term in terms {
				result += term.coefficient*pow(x, Double(term.degree))
			}
			return result
		}
		public func makeDerivative() -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree > 0 {
					newTerms.append(Term(term.coefficient*Double(term.degree), exp: term.degree-1))
				}
			}
			if newTerms.isEmpty { return PolynomialEQ.zero }
			return PolynomialEQ(terms: newTerms)
		}
		public func integrate(plus c: Double) -> MathEquation {
			var newTerms = [Term]()
			for term in terms {
				if term.degree >= 0 {
					newTerms.append(Term(term.coefficient/Double((term.degree+1)), exp: term.degree+1))
				}
			}
			newTerms.append(Term(c, exp: 0))
			return PolynomialEQ(terms: newTerms)
		}
		public mutating func organizeTerms() {
			let highestDegree = terms.map {$0.degree}.max() ?? 0
			let lowestDegree = terms.map {$0.degree}.min() ?? 0
			var tempTerms = [Term]()
			for i in lowestDegree...highestDegree {
				let filteredTerms = terms.filter({$0.degree == i})
				if filteredTerms.count == 0 { continue }
				if filteredTerms.count == 1 {
					tempTerms.append(filteredTerms.first!)
				} else {
					if let term = filteredTerms.merge(combining: { t1, t2 in
						Term(t1.coefficient+t2.coefficient, exp: i)
					}) {
						tempTerms.append(term)
					}
				}
			}
			self.terms = tempTerms.sorted(by: {$0.degree > $1.degree})
		}
		func makeMonic() -> PolynomialEQ {
			// Divides the polynomial with the leading coefficient to make the polynomial monic
			guard terms.count > 0 else { return self }
			var monicTerms = terms
			for i in terms.indices {
				monicTerms[i] = terms[i]/(terms[0].coefficient)
			}
			return PolynomialEQ(terms: monicTerms)
		}
		
		public mutating func scale(by rhs: Double) {
			for i in terms.indices {
				terms[i] = terms[i]*rhs
			}
		}
		public func scaled(by rhs: Double) -> Self {
			var terms: [Term] = []
			for t in self.terms {
				terms.append(t*rhs)
			}
			return Self.init(terms: terms)
		}
		public func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double = 120) -> (min: Double, max: Double) {
			var minValue, maxValue: Double
			if let range {
				minValue = min(callAsFunction(range.lowerBound),callAsFunction(range.upperBound))
				maxValue = max(callAsFunction(range.lowerBound),callAsFunction(range.upperBound))
			} else {
				if isEvenFunction {
					if terms.max(by: {$0.degree < $1.degree})!.coefficient < 0 {
						minValue = -Double.infinity
						maxValue = callAsFunction(0)
					} else {
						minValue = callAsFunction(0)
						maxValue = Double.infinity
					}
				} else {
					return (-Double.infinity, Double.infinity)
				}
			}
			for x in criticalPoints(in: range) {
				minValue = min(minValue, callAsFunction(x))
				maxValue = max(maxValue, callAsFunction(x))
			}
			return (minValue, maxValue)
		}
		public func criticalPoints(in range: ClosedRange<Double>?) -> [Double] {
			let dydx = makeDerivative() as! PolynomialEQ
			return (try? dydx.zeros(in: range)) ?? []
		}
		public static let zero: Math.PolynomialEQ = Math.PolynomialEQ(y: 0)
		public static func + (lhs: Math.PolynomialEQ, rhs: Math.PolynomialEQ) -> Math.PolynomialEQ {
			Math.PolynomialEQ(terms: lhs.terms+rhs.terms)
		}
		public static func - (lhs: Math.PolynomialEQ, rhs: Math.PolynomialEQ) -> Math.PolynomialEQ {
			let rhsTerms = rhs.terms.map({Term(-$0.coefficient, exp: $0.degree)})
			return Math.PolynomialEQ(terms: lhs.terms+rhsTerms)
		}
		
		public var description: String {
			let termDescriptions = terms.map { term in
				term.description
			}
			return termDescriptions.joined(separator: " + ")
		}
		public struct Term: CustomStringConvertible, Sendable, Equatable {

			private(set) var degree: UInt
			private(set) var coefficient: Double
			
			public init(_ coefficient: Double, exp degree: UInt = 0) {
				self.degree = degree
				self.coefficient = coefficient
			}
			public static func + (_ lhs: Self, _ rhs: Self) -> Self {
				return Term(lhs.coefficient+rhs.coefficient, exp: lhs.degree)
			}
			public static func * (_ lhs: Self, _ rhs: Double) -> Self {
				return Term(lhs.coefficient*rhs, exp: lhs.degree)
			}
			public static func / (_ lhs: Self, _ rhs: Double) -> Self {
				return Term(lhs.coefficient/rhs, exp: lhs.degree)
			}
			public static func * (_ lhs: Double, _ rhs: Self) -> Self {
				return Term(lhs*rhs.coefficient, exp: rhs.degree)
			}
			
			public var description: String {
				if degree == 1 {
					return "\(coefficient)x"
				} else if degree == 0 {
					return "\(coefficient)"
				} else {
					return "\(coefficient)x^\(degree)"
				}
			}
		}
		
		public var isEvenFunction: Bool {
			guard !terms.isEmpty else { return true }
			let d = terms.map({$0.degree}).max()!
			return d % 2 == 0
		}
		public func zeros(in range: ClosedRange<Double>?) throws -> [Double] {
			guard !terms.isEmpty else { throw PolynomialEQError.emptyPolynomial }
			let degree = terms.map({$0.degree}).max() ?? 0

			// Handle low degree polynomials
			if degree == 0  {
				if callAsFunction(0) == 0 {
					if range?.contains(0) ?? true { return [0] } else { return [] }
				} else { return [] }
			}
			if degree == 1 {
				let a = terms.first(where: { $0.degree == 1 })?.coefficient ?? 0
				let b = terms.first(where: { $0.degree == 0 })?.coefficient ?? 0
				if a == 0 && b == 0 {
					if range?.contains(0) ?? true { return [0] } else { return [] }
				} else if a == 0 { return [] }
				if range?.contains(-b/a) ?? true { return [-b/a] } else { return [] }
			}
			var zeros = [Double]()
			if degree == 2 {
				let a = terms.first(where: { $0.degree == 2 })?.coefficient ?? 0
				let b = terms.first(where: { $0.degree == 1 })?.coefficient ?? 0
				let c = terms.first(where: { $0.degree == 0 })?.coefficient ?? 0
				let discriminant = b*b-4*a*c
				guard a != 0, discriminant >= 0 else { return [] }
				if discriminant == 0 {
					zeros = [(-b)/(2 * a)]
				} else {
					let sqrtDisc = sqrt(discriminant)
					zeros = [(-b + sqrtDisc) / (2 * a), (-b - sqrtDisc) / (2 * a)]
				}
			}
			// Handle higher order polynomials
			let companion = try companionMatrix()
			let eigenvalues = try companion.eigenvalues()
			
			// Filter real roots
			zeros = eigenvalues.compactMap({
				if $0.i.isApproxEqual(to: 0) {
					return $0.x
				} else {
					return nil
				}
			})
			// Verify zeros
			for zero in zeros {
				if !callAsFunction(zero).isApproxEqual(to: 0) { throw PolynomialEQError.invalidZeros }
			}
			// Filter zeros based on range
			if let range {
				zeros = zeros.filter({range.contains($0)})
			}
			return zeros
		}
		private func polynomialCoefficients() -> [Double] {
			let maxDegree = terms.map(\.degree).max() ?? 0
			var coeffs = Array(repeating: 0.0, count: Int(maxDegree) + 1)
			for term in terms {
				coeffs[Int(term.degree)] += term.coefficient
			}
			return coeffs.reversed() // Highest degree first
		}
		private func companionMatrix() throws -> Matrix<Double> {
			let coeffs = polynomialCoefficients()
			let n = coeffs.count - 1
			guard n > 0 else { throw PolynomialEQError.emptyPolynomial }
			guard coeffs[0] != 0 else { throw PolynomialEQError.notMonic }
			
			// Normalize to reversed monic polynomial
			let a = coeffs.map({ $0 / coeffs[0] }).reversedRawCollection()
			
			// Build companion matrix
			var companion = Matrix<Double>(size: n)
			for i in 1..<n {
				companion[i,i-1] = 1.0
			}
			for i in 0..<n {
				companion[i,n-1] = -a[i]
			}
			return companion
		}
	}
	
	public enum PolynomialEQError: Error {
		case notMonic, emptyPolynomial, invalidZeros
	}
	
	/// A defined value at an unspecified location.
	public struct Impulse: MathEquation {
		
		public private(set) var term: Double

		public init(term: Double) {
			self.term = term
		}
		public func callAsFunction(_ x: Double) -> Double {
			return term
		}
		public func makeDerivative() -> MathEquation {
			return NullEQ()
		}
		public func integrate(plus c: Double) -> MathEquation {
			return PolynomialEQ(terms: PolynomialEQ.Term(term, exp: 0))
		}
		public mutating func scale(by rhs: Double) {
			term *= rhs
		}
		public func scaled(by rhs: Double) -> Self {
			return Impulse(term: term*rhs)
		}
		public func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double = 120) -> (min: Double, max: Double) {
			let minValue = min(0,term)
			let maxValue = max(0,term)
			return (minValue, maxValue)
		}
		public func criticalPoints(in: ClosedRange<Double>?) -> [Double] {
			return []
		}
		public var description: String {
			return "\(term) (impulse)"
		}
	}
	
	/// Equivalent to zero
	public struct NullEQ: MathEquation {
		
		public init() {}
		public func callAsFunction(_ x: Double) -> Double {
			return 0
		}
		public func makeDerivative() -> MathEquation {
			return NullEQ()
		}
		public func integrate(plus c: Double) -> MathEquation {
			return PolynomialEQ(terms: PolynomialEQ.Term(c, exp: 0))
		}
		public mutating func scale(by rhs: Double) {}
		public func scaled(by rhs: Double) -> Self {
			return self
		}
		public var description: String {
			return "0 (null)"
		}
		public func minmaxValues(in range: ClosedRange<Double>?, num_dx: Double = 120) -> (min: Double, max: Double) {
			return (0, 0)
		}
		public func criticalPoints(in: ClosedRange<Double>?) -> [Double] {
			return []
		}
	}
	
}
//public struct AnyMathEquation: MathEquation {
//	private let _callAsFunction: (Double) -> Double
//	private let _makeDerivative: () -> MathEquation
//	private let _integrate: (Double) -> MathEquation
//	private let _description: () -> String
//
//	public init<E: MathEquation>(_ equation: E) {
//		_callAsFunction = equation.callAsFunction
//		_makeDerivative = equation.makeDerivative
//		_integrate = equation.integrate
//		_description = { equation.description }
//	}
//
//	public func callAsFunction(_ x: Double) -> Double {
//		_callAsFunction(x)
//	}
//
//	public func makeDerivative() -> MathEquation {
//		_makeDerivative()
//	}
//
//	public func integrate(plus c: Double) -> MathEquation {
//		_integrate(c)
//	}
//
//	public var description: String {
//		_description()
//	}
//}
