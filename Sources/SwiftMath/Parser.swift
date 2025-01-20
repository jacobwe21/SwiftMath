//
//  Parser.swift
//  MySwiftMath
//
//  Created by Jacob W Esselstyn on 12/30/24.
// 	Based on github.com/peredaniel/MathExpression. Copyright © 2019 Pedro Daniel Prieto Martínez. Distributed under MIT License.
//

import Foundation

public struct MathExpression {
	public enum ValidationError: Error {
		case emptyExpression
		case misplacedBrackets
		case unevenOpeningClosingBracketNumber
		case invalidConsecutiveOperators(String)
		case startsWithNonSumOrSubtractionOperator(String)
		case endsWithOperator(String)
		case tooManyDecimals
		case notMath
	}

	private let formula: MathFormula

	public init(
		_ string: String,
		transformation: @escaping (String) -> Double = { Double($0) ?? .zero }
	) throws {
		formula = try MathFormula(string, transformation: transformation)
	}

	init(_ formula: MathFormula) {
		self.formula = formula
	}

	public func evaluate() -> Double {
		switch formula.evaluationState() {
		case .isNumeric(let value):
			return value
		case .startsWithSymbol(let symbol):
			switch symbol {
			case .sum:
				return MathExpression(formula.dropingInitialValue()).evaluate()
			case .subtraction:
				return MathExpression(formula.replaceSubtractionByNegative()).evaluate()
			}
		case .containsBracket(let brackets):
			return MathExpression(formula.evaluatingExpression(between: brackets)).evaluate()
		case .canApplyOperator(let mathOperator):
			return mathOperator.apply(
				to: formula.decompose(with: mathOperator).map { MathExpression($0) }.evaluate()
			)
		case .canApplyTransformation:
			return formula.applyTransformation()
		}
	}
}

fileprivate extension Array where Element == MathExpression {
	func evaluate() -> [Double] {
		map { $0.evaluate() }
	}
}

struct MathFormula {
	enum EvaluationState {
		case isNumeric(Double)
		case startsWithSymbol(AdditiveOperator)
		case containsBracket(MathBrackets)
		case canApplyOperator(MathOperator)
		case canApplyTransformation
	}

	private let string: String
	private let transformation: (String) -> Double

	init(
		_ string: String,
		transformation: @escaping (String) -> Double = { Double($0) ?? .zero }
	) throws {
		guard string.containsOnlyMath else {
			let notDigitsOrSymbols = NSCharacterSet.decimalDigits.union(CharacterSet(charactersIn: ".,+/xX*--֊־‐‑‒–—―－eE(){}[]")).union(NSCharacterSet.whitespaces).inverted
			let badCharater: Character? = string.first { !notDigitsOrSymbols.contains($0.unicodeScalars.first!) }
			print(badCharater ?? "No bad characters")
			throw MathExpression.ValidationError.notMath
		}
		self.string = string.cleanMathExpression()
		var decimalCount: Int = 0
		for i in self.string {
			if i == "." { decimalCount += 1 }
			if decimalCount > 1 { throw MathExpression.ValidationError.tooManyDecimals }
		}
		self.transformation = transformation
		try validate()
	}

	init(
		validString string: String,
		transformation: @escaping (String) -> Double
	) {
		self.string = string.polishMathExpression()
		self.transformation = transformation
	}

	func evaluationState(validating: Bool = false) -> EvaluationState {
		if let value = Double(string) {
			return .isNumeric(value)
		}

		for additiveOperator in AdditiveOperator.allCases {
			if starts(with: additiveOperator) {
				return .startsWithSymbol(additiveOperator)
			}
		}

		for brackets in MathBrackets.allCases {
			if containsBracket(brackets) {
				return .containsBracket(brackets)
			}
		}

		if let firstOperator = getFirstOperator(validating: validating) {
			return .canApplyOperator(firstOperator)
		}

		return .canApplyTransformation
	}

	fileprivate func getFirstOperator(validating: Bool = false) -> MathOperator? {
		let cases = validating ? MathOperator.validationCases : MathOperator.evaluationCases
		return cases.first { string.contains($0.rawValue) }
	}
	
	func applyTransformation() -> Double {
		transformation(string)
	}

	fileprivate func decompose(with mathOperator: MathOperator) -> [MathFormula] {
		var finalString = string
		if let _ = MathOperator.validConsecutiveOperatorsDuringEvaluation.first(where: { string.contains($0.key) }) {
			for (doubleOperator, combinedOperator) in MathOperator.validConsecutiveOperatorsDuringEvaluation {
				finalString = finalString.replacingOccurrences(of: doubleOperator, with: combinedOperator)
			}
		}
		return finalString.split(separator: mathOperator.character).map {
			MathFormula(validString: String($0), transformation: transformation)
		}
	}

	fileprivate func evaluatingExpression(between bracket: MathBrackets) -> MathFormula {
		let stringWithBrackets = (try? getString(between: bracket)) ?? ""

		let value = MathExpression(
			MathFormula(
				validString: String(stringWithBrackets.dropFirst().dropLast()),
				transformation: transformation
			)
		).evaluate()
		return MathFormula(
			validString: string.replacingOccurrences(of: stringWithBrackets, with: value.avoidScientificNotation()),
			transformation: transformation
		)
	}

	fileprivate func dropingInitialValue() -> MathFormula {
		MathFormula(
			validString: String(string.dropFirst()),
			transformation: transformation
		)
	}
	fileprivate func replaceSubtractionByNegative() -> MathFormula {
		return MathFormula(
			validString: MathOperator.negative.rawValue + String(string.dropFirst()),
			transformation: transformation
		)
	}
}

fileprivate extension MathFormula {
	func containsBracket(_ bracket: MathBrackets) -> Bool {
		string.contains(bracket.opening) || string.contains(bracket.closing)
	}

	func getString(between bracket: MathBrackets) throws -> String {
		try String(string.map { $0 }.characters(between: bracket))
	}

	func starts(with additiveOperator: AdditiveOperator) -> Bool {
		string.first == additiveOperator.character
	}

	func validate() throws {
		guard !string.isEmpty else {
			throw MathExpression.ValidationError.emptyExpression
		}

		try string.validateStartingAndEndingCharacters()
		try string.validateNoInvalidConsecutiveOperators()

		switch evaluationState(validating: true) {
		case .isNumeric:
			break
		case .startsWithSymbol(let symbol):
			switch symbol {
			case .sum:
				_ = try MathFormula(String(string.dropFirst()), transformation: transformation)
			case .subtraction:
				_ = try MathFormula(MathOperator.negative.rawValue + String(string.dropFirst()), transformation: transformation)
			}
		case .containsBracket(let brackets):
			let stringWithBrackets = try getString(between: brackets)
			_ = try MathFormula(
				String(stringWithBrackets.dropFirst().dropLast()),
				transformation: transformation
			)
			_ = try MathFormula(
				string.replacingOccurrences(of: stringWithBrackets, with: "0"),
				transformation: transformation
			)

		case .canApplyOperator(let mathOperator):
			_ = try string.split(separator: mathOperator.character).map {
				try MathFormula(String($0), transformation: transformation)
			}
		case .canApplyTransformation:
			break
		}
	}
}

// MARK: - Helper extensions

fileprivate extension String {
	func cleanMathExpression() -> String {
		var str = self
		str = str.replacingOccurrences(of: "֊", with: "-")
		str = str.replacingOccurrences(of: "־", with: "-")
		str = str.replacingOccurrences(of: "᠆", with: "-")
		str = str.replacingOccurrences(of: "‐", with: "-")
		str = str.replacingOccurrences(of: "‑", with: "-")
		str = str.replacingOccurrences(of: "‒", with: "-")
		str = str.replacingOccurrences(of: "–", with: "-")
		str = str.replacingOccurrences(of: "—", with: "-")
		str = str.replacingOccurrences(of: "―", with: "-")
		str = str.replacingOccurrences(of: "－", with: "-")
		str = str.replacingOccurrences(of: "X", with: "*")
		str = str.replacingOccurrences(of: "x", with: "*")
		str = str.replacingOccurrences(of: "{", with: "(")
		str = str.replacingOccurrences(of: "}", with: ")")
		str = str.replacingOccurrences(of: "[", with: "(")
		str = str.replacingOccurrences(of: "]", with: ")")
		str = str.replacingOccurrences(of: ")(", with: ")*(")
		if Locale.current.decimalSeparator == "." {
			str = str.replacingOccurrences(of: ",", with: "")
		} else {
			str = str.replacingOccurrences(of: ".", with: "")
			str = str.replacingOccurrences(of: ",", with: ".")
		}
		return str.polishMathExpression()
	}
	func polishMathExpression() -> String {
		var finalString = self.removingWhitespace //replacingOccurrences(of: " ", with: "")
		for (doubleOperator, combinedOperator) in MathOperator.validConsecutiveOperators {
			finalString = finalString.replacingOccurrences(of: doubleOperator, with: combinedOperator)
		}
		return finalString
	}

	func validateNoInvalidConsecutiveOperators() throws {
		for consecutiveOperators in MathOperator.invalidConsecutiveOperators {
			if let _ = range(of: consecutiveOperators) {
				throw MathExpression.ValidationError.invalidConsecutiveOperators(consecutiveOperators)
			}
		}
	}

	func validateStartingAndEndingCharacters() throws {
		guard let first = first, let last = last else { return }
		if MathOperator.multiplicativeOperators.map({ $0.rawValue }).contains(String(first)) {
			throw MathExpression.ValidationError.startsWithNonSumOrSubtractionOperator(String(first))
		}

		if MathOperator.validationCases.map({ $0.rawValue }).contains(String(last)) {
			throw MathExpression.ValidationError.endsWithOperator(String(last))
		}
	}
}

fileprivate extension Array where Element == String.Element {
	func characters(between brackets: MathBrackets) throws -> [Character] {
		guard let initialIndex = lastIndex(of: Character(brackets.opening)) else {
			throw MathExpression.ValidationError.unevenOpeningClosingBracketNumber
		}

		guard let finalIndex = self[initialIndex..<count].firstIndex(of: Character(brackets.closing)) else {
			throw MathExpression.ValidationError.misplacedBrackets
		}

		return Array(self[initialIndex...finalIndex])
	}
}

extension MathExpression.ValidationError: Equatable {
	public static func == (lhs: MathExpression.ValidationError, rhs: MathExpression.ValidationError) -> Bool {
		switch (lhs, rhs) {
		case (.emptyExpression, .emptyExpression),
			(.misplacedBrackets, .misplacedBrackets),
			(.unevenOpeningClosingBracketNumber, .unevenOpeningClosingBracketNumber):
			return true
		case (.invalidConsecutiveOperators(let lhsValue), .invalidConsecutiveOperators(let rhsValue)),
			(.startsWithNonSumOrSubtractionOperator(let lhsValue), .startsWithNonSumOrSubtractionOperator(let rhsValue)),
			(.endsWithOperator(let lhsValue), .endsWithOperator(let rhsValue)):
			return lhsValue == rhsValue
		default:
			return false
		}
	}
}


// MARK: Parser Model Structure
enum MathBrackets: CaseIterable {
	case parenthesis//, bracket, curly
	var opening: String {
		switch self {
		case .parenthesis: 	return "("
		//case .bracket: 		return "["
		//case .curly: 		return "{"
		}
	}
	var closing: String {
		switch self {
		case .parenthesis: 	return ")"
		//case .bracket:		return "]"
		//case .curly:		return "}"
		}
	}
}
enum AdditiveOperator: String, CaseIterable {
	case sum = "+"
	case subtraction = "-"

	var character: Character {
		Character(rawValue)
	}
}
enum MathOperator: String {
	case negative = "_"
	case product = "*"
	case division = "/"
	case sum = "+"
	case subtraction = "-"
	
	var character: Character {
		Character(rawValue)
	}

	func apply(to args: [Double]) -> Double {
		switch self {
		case .sum: return args.sum()
		case .subtraction: return args.subtract()
		case .product: return args.multiply()
		case .division: return args.divide()
		case .negative: return args.first?.negative ?? .zero
		}
	}

	static var evaluationCases: [MathOperator] {
		[.sum, .subtraction, .product, .division, .negative]
	}

	static var multiplicativeOperators: [MathOperator] {
		[.product, .division]
	}

	static var validationCases: [MathOperator] {
		[.product, .division, .sum, .subtraction]
	}
	
	static var validConsecutiveOperators: [String: String] {
		[
			MathOperator.sum.rawValue + MathOperator.subtraction.rawValue: MathOperator.subtraction.rawValue,
			MathOperator.subtraction.rawValue + MathOperator.sum.rawValue: MathOperator.subtraction.rawValue,
			MathOperator.sum.rawValue + MathOperator.sum.rawValue: MathOperator.sum.rawValue,
			MathOperator.subtraction.rawValue + MathOperator.subtraction.rawValue: MathOperator.sum.rawValue
		]
	}

	static var validConsecutiveOperatorsDuringEvaluation: [String: String] {
		[
			MathOperator.negative.rawValue + MathOperator.sum.rawValue: MathOperator.negative.rawValue,
			MathOperator.negative.rawValue + MathOperator.subtraction.rawValue: MathOperator.sum.rawValue,
			MathOperator.sum.rawValue + MathOperator.subtraction.rawValue: MathOperator.sum.rawValue + MathOperator.negative.rawValue,
			MathOperator.product.rawValue + MathOperator.subtraction.rawValue: MathOperator.product.rawValue + MathOperator.negative.rawValue,
			MathOperator.division.rawValue + MathOperator.subtraction.rawValue: MathOperator.division.rawValue + MathOperator.negative.rawValue
		]
	}

	static var invalidConsecutiveOperators: [String] {
		[
			MathOperator.sum.rawValue + MathOperator.product.rawValue,
			MathOperator.sum.rawValue + MathOperator.division.rawValue,
			MathOperator.subtraction.rawValue + MathOperator.product.rawValue,
			MathOperator.subtraction.rawValue + MathOperator.division.rawValue,
			MathOperator.product.rawValue + MathOperator.product.rawValue,
			MathOperator.product.rawValue + MathOperator.division.rawValue,
			MathOperator.division.rawValue + MathOperator.product.rawValue,
			MathOperator.division.rawValue + MathOperator.division.rawValue
		]
	}
}

// MARK: - Helper extensions
fileprivate extension Array where Element == Double {
	func sum() -> Element {
		guard let last = last else { return .zero }
		return last + dropLast().sum()
	}

	func subtract() -> Double {
		guard let first = first else { return .zero }
		return first + reversed().dropLast().map { $0.negative }.sum()
	}

	func multiply() -> Double {
		guard !contains(.zero) else { return .zero }
		guard let last = last else { return 1.0 }
		return last * dropLast().multiply()
	}

	func divide() -> Double {
		guard let numerator = first else { return 1.0 }
		return numerator / reversed().dropLast().multiply()
	}
}
fileprivate extension FloatingPoint {
	var negative: Self {
		var negativeValue = self
		negativeValue.negate()
		return negativeValue
	}

	func avoidScientificNotation() -> String {
		Formatter.avoidScientificNotation.string(for: self) ?? ""
	}
}
fileprivate extension Formatter {
	static let avoidScientificNotation: NumberFormatter = {
		let numberFormatter = NumberFormatter()
		numberFormatter.maximumFractionDigits = 16
		numberFormatter.numberStyle = .decimal
		numberFormatter.decimalSeparator = "."
		numberFormatter.usesGroupingSeparator = false
		return numberFormatter
	}()
}
extension String {
	var containsOnlyMath: Bool {
		let notDigitsOrSymbols = NSCharacterSet.decimalDigits.union(CharacterSet(charactersIn: "+--֊־‐‑‒–—―－/xX*.,(){}[]eE_")).union(NSCharacterSet.whitespaces).inverted
		return rangeOfCharacter(from: notDigitsOrSymbols, options: String.CompareOptions.literal, range: nil) == nil
	}
}

