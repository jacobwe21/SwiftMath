//
//  Measurements.swift
//
//  Created by Jacob W Esselstyn on 6/23/23.
//

import Foundation
import MySwift

public protocol EngineeringUnit: Dimension {
	associatedtype EngDimension: EngineeringUnit
	static var allEngineeringUnits: [EngDimension] { get }
	static var allEngineeringUnitSymbols: [String] { get }
	static var allImperialEngineeringUnitSymbols: [String] { get }
	static var allSIEngineeringUnitSymbols: [String] { get }
	//static var deafaultImperialUnit: EngDimension { get }
	//static var deafaultSIUnit: EngDimension { get }
	var positiveOnly: Bool { get }
	var isImperial: Bool { get }
}
public extension EngineeringUnit {
	var positiveOnly: Bool { false }
	static var allEngineeringUnitSymbols: [String] { allEngineeringUnits.map({$0.symbol}) }
	static var allImperialEngineeringUnitSymbols: [String] { allEngineeringUnits.filter({$0.isImperial}).map({$0.symbol}) }
	static var allSIEngineeringUnitSymbols: [String] { allEngineeringUnits.filter({!$0.isImperial}).map({$0.symbol}) }
	//static var deafaultSIUnit: EngDimension { baseUnit() }
}

///  The inverse (1/x) of `UnitTemperature`
public class UnitInverseTemperature: Dimension, EngineeringUnit {
	public static let allEngineeringUnits: [UnitInverseTemperature] = [.inverseKelvin, .inverseCelsius, .inverseFahrenheit]
	
	public static let inverseKelvin = UnitInverseTemperature(symbol: "1/K", converter: UnitConverterLinear(coefficient: 1))
	public static let inverseCelsius = UnitInverseTemperature(symbol: "1/°C", converter: UnitConverterInverting(coefficient: 1, constant: -237.15))
	public static let inverseFahrenheit = UnitInverseTemperature(symbol: "1/°F", converter: UnitConverterInverting(coefficient: 5/9, constant: -457.87))

	public override class func baseUnit() -> Self {
		return UnitInverseTemperature.inverseKelvin as! Self
	}
	
	override public static var supportsSecureCoding: Bool { true }
	public var isImperial: Bool {
		if self ==|| [.inverseFahrenheit] {
			return true
		} else { return false }
	}
	public var positiveOnly: Bool { self == .inverseKelvin }
	
	class UnitConverterInverting: UnitConverter {
		var coefficient: Double
		var constant: Double
		var doubleInversion: Bool
		init(doubleInversion: Bool = true, coefficient: Double = 1, constant: Double = 0) {
			self.coefficient = coefficient
			self.constant = constant
			self.doubleInversion = doubleInversion
		}
		override func baseUnitValue(fromValue value: Double) -> Double {
			let x: Double
			if doubleInversion {
				guard value != 0 else { return .infinity }
				x = 1/value
			} else { x = value }
			let x2 = (x - constant)
			guard x2 != 0 else { return .infinity }
			return coefficient/x2
		}
		override func value(fromBaseUnitValue baseUnitValue: Double) -> Double {
			guard baseUnitValue != 0 else { return .infinity }
			let x = coefficient/baseUnitValue + constant
			if doubleInversion {
				guard x != 0 else { return .infinity }
				return 1/x
			} else {
				return x
			}
		}
	}
}
///  A unit of measure for density (technically the same as `UnitConcentrationMass`, but kg/m³ is the base unit for density).
public class UnitDensity: Dimension, EngineeringUnit {

	public static let allEngineeringUnits: [UnitDensity] = [.kilogramPerCubicMeter, .poundsPerCubicFoot]
	
	public static let kilogramPerCubicMeter = UnitDensity(symbol: "kg/m³", converter: UnitConverterLinear(coefficient: 1.0))
	
	public static let poundsPerCubicFoot = UnitDensity(symbol: "lb/ft³", converter: UnitConverterLinear(coefficient: 16.018463))
	
	public override class func baseUnit() -> Self {
		return UnitDensity.kilogramPerCubicMeter as! Self
	}
	override public static var supportsSecureCoding: Bool { true }
	public var positiveOnly: Bool { true }
	
	public var isImperial: Bool {
		if self ==|| [.poundsPerCubicFoot] {
			return true
		} else { return false }
	}
}
///  A unit of measure for force
public class UnitForce: Dimension, EngineeringUnit {
	public static let newton = UnitForce(symbol: "N", converter: UnitConverterLinear(coefficient: 1.0))
	public static let kilonewton = UnitForce(symbol: "kN", converter: UnitConverterLinear(coefficient: 1000))
	public static let pound = UnitForce(symbol: "lb", converter: UnitConverterLinear(coefficient: 4.44822162825))
	public static let kip = UnitForce(symbol: "k", converter: UnitConverterLinear(coefficient: 4448.22162825))
	
	public override class func baseUnit() -> Self {
		return UnitForce.newton as! Self
	}
	public static let allEngineeringUnits: [UnitForce] = [.newton,.pound,.kip,.kilonewton]
	override public static var supportsSecureCoding: Bool { true }
	public var isImperial: Bool {
		if self ==|| [.pound,.kip] {
			return true
		} else { return false }
	}
}
///  A unit of measure for force distributed over a linear area
public class UnitLinearForce: Dimension, EngineeringUnit {
	public static let newtonsPerMeter = UnitLinearForce(symbol: "N/m", converter: UnitConverterLinear(coefficient: 1.0))
	public static let kilonewtonsPerMeter = UnitLinearForce(symbol: "kN/m", converter: UnitConverterLinear(coefficient: 1000))
	public static let newtonsPerCentimeter = UnitLinearForce(symbol: "N/cm", converter: UnitConverterLinear(coefficient: 100.0))
	public static let newtonsPerMillimeter = UnitLinearForce(symbol: "N/mm", converter: UnitConverterLinear(coefficient: 1000.0))
	public static let kilonewtonsPerMillimeter = UnitLinearForce(symbol: "kN/mm", converter: UnitConverterLinear(coefficient: 1_000_000.0))
	
	public static let poundsPerFoot = UnitLinearForce(symbol: "lb/ft", converter: UnitConverterLinear(coefficient: 14.5939029))
	public static let kipsPerFoot = UnitLinearForce(symbol: "k/ft", converter: UnitConverterLinear(coefficient: 14593.9029))
	public static let poundsPerInch = UnitLinearForce(symbol: "lb/in", converter: UnitConverterLinear(coefficient: 175.126835))
	public static let kipsPerInch = UnitLinearForce(symbol: "k/in", converter: UnitConverterLinear(coefficient: 175126.835))
	
	public override class func baseUnit() -> Self {
		return UnitLinearForce.newtonsPerMeter as! Self
	}
	public static let allEngineeringUnits: [UnitLinearForce] = [.newtonsPerMeter,.kilonewtonsPerMeter,.newtonsPerCentimeter,.newtonsPerMillimeter,.kilonewtonsPerMillimeter,.poundsPerFoot,.poundsPerInch,.kipsPerFoot,.kipsPerInch]
	override public static var supportsSecureCoding: Bool { true }
	public var isImperial: Bool {
		if self ==|| [.poundsPerFoot,.poundsPerInch,.kipsPerInch,.kipsPerFoot] {
			return true
		} else { return false }
	}
}
/// A unit of measure for work (force times distance) - Technically the same as the UnitEnergy class
public class UnitWork: Dimension, EngineeringUnit {
	public static let newtonMeters = UnitWork(symbol: "N•m", converter: UnitConverterLinear(coefficient: 1.0))
	public static let kilonewtonMeters = UnitWork(symbol: "kN•m", converter: UnitConverterLinear(coefficient: 1000.0))
	public static let kilonewtonMillimeters = UnitWork(symbol: "kN•mm", converter: UnitConverterLinear(coefficient: 1.0))
	public static let newtonCentimeters = UnitWork(symbol: "N•cm", converter: UnitConverterLinear(coefficient: 0.01))
	public static let newtonMillimeters = UnitWork(symbol: "N•mm", converter: UnitConverterLinear(coefficient: 0.001))
	
	public static let poundFeet = UnitWork(symbol: "lb•ft", converter: UnitConverterLinear(coefficient: 1.35582))
	public static let kipFeet = UnitWork(symbol: "k•ft", converter: UnitConverterLinear(coefficient: 1355.8179483314))
	public static let kipInches = UnitWork(symbol: "k•in", converter: UnitConverterLinear(coefficient: 112.9848293))
	public static let poundInches = UnitWork(symbol: "lb•in", converter: UnitConverterLinear(coefficient: 0.11298482933))
	
	public override class func baseUnit() -> Self {
		return UnitWork.newtonMeters as! Self
	}
	public static let allEngineeringUnits: [UnitWork] = [ .newtonMeters, .newtonCentimeters, .newtonMillimeters, .kilonewtonMeters, .kilonewtonMillimeters, .poundFeet, .kipFeet, .kipInches, .poundInches]
	override public static var supportsSecureCoding: Bool { true }
	public var isImperial: Bool {
		if self ==|| [.poundFeet,.poundInches,.kipInches,.kipFeet] {
			return true
		} else { return false }
	}
}
extension UnitLength: EngineeringUnit {
	public static let allEngineeringUnits: [UnitLength] = [.meters,.centimeters,.millimeters,.feet,.inches]
//	public static let allEngineeringUnits: [UnitLength] = [.kilometers,.meters,.centimeters,.millimeters,.miles,.yards,.feet,.inches]
	
	public var isImperial: Bool {
		if self ==|| [.inches,.feet,.yards,.miles,.furlongs,.fathoms] {
			return true
		} else { return false }
	}
}
public extension Measurement where UnitType: UnitLength {
	func squared() -> Measurement<UnitArea> {
		Self.multiplyMeasurementsIntoNewMeasurement(a: self, b: self, returnType: UnitArea.self)
	}
	func cubed() -> Measurement<UnitVolume> {
		Self.multiplyMeasurementsIntoNewMeasurement(a: self.squared(), b: self, returnType: UnitVolume.self)
	}
	func tesseracted() -> Measurement<UnitTesseract> {
		Self.multiplyMeasurementsIntoNewMeasurement(a: self.squared(), b: self.squared(), returnType: UnitTesseract.self)
	}
}

extension UnitArea: EngineeringUnit {
	public static let allEngineeringUnits: [UnitArea] = [.squareMeters,.squareCentimeters,.squareMillimeters,.squareFeet,.squareInches]
//	public static let allEngineeringUnits: [UnitArea] = [.squareKilometers,.squareMeters,.squareCentimeters,.squareMillimeters,.squareMiles,.squareYards,.squareFeet,.squareInches]
	
	/// The unit that results by taking the square root of a UnitArea unit.
	var linearBaseUnit: UnitLength? {
		switch self {
		case .squareMegameters: UnitLength.megameters
		case .squareKilometers: UnitLength.kilometers
		case .squareMeters: UnitLength.meters
		case .squareCentimeters: UnitLength.centimeters
		case .squareMillimeters: UnitLength.millimeters
		case .squareMicrometers: UnitLength.micrometers
		case .squareNanometers: UnitLength.nanometers
		case .squareMiles: UnitLength.miles
		case .squareYards: UnitLength.yards
		case .squareFeet: UnitLength.feet
		case .squareInches: UnitLength.inches
		default: nil
		}
	}
	public var isImperial: Bool {
		if self ==|| [.acres,.squareFeet,.squareInches,.squareMiles,.squareYards] {
			return true
		} else { return false }
	}
}
public extension Measurement where UnitType: UnitArea {
	func sqrt() -> Measurement<UnitLength> {
		if let linearUnit = self.unit.linearBaseUnit {
			return Measurement<UnitLength>(value: Darwin.sqrt(self.value), unit: linearUnit)
		} else {
			return Measurement<UnitLength>(value: Darwin.sqrt(self.baseUnitValue()), unit: .baseUnit())
		}
	}
	func squared() -> Measurement<UnitTesseract> {
		Self.multiplyMeasurementsIntoNewMeasurement(a: self, b: self, returnType: UnitTesseract.self)
	}
}

extension UnitVolume: EngineeringUnit {
	public static let allEngineeringUnits: [UnitVolume] = [.cubicMeters,.cubicCentimeters,.cubicMillimeters,.cubicFeet,.cubicInches]
//	public static let allEngineeringUnits: [UnitVolume] = [.cubicMeters,.cubicCentimeters,.cubicMillimeters,.cubicYards,.cubicFeet,.cubicInches,.gallons,.liters]
	public var isImperial: Bool {
		if self ==|| [.acreFeet,.bushels,.cubicFeet,.cubicInches,.cubicMiles,.cubicYards,.cups,.imperialFluidOunces,.imperialPints,.imperialQuarts,.imperialGallons,.imperialTeaspoons,.imperialTablespoons] {
			return true
		} else { return false }
	}
}
public class UnitTesseract: Dimension, EngineeringUnit {
	public static let tesseractMeters = UnitTesseract(symbol: "m⁴", converter: UnitConverterLinear(coefficient: 1))
	public static let tesseractCentimeters = UnitTesseract(symbol: "cm⁴", converter: UnitConverterLinear(coefficient: 0.00000001))
	public static let tesseractMillimeters = UnitTesseract(symbol: "mm⁴", converter: UnitConverterLinear(coefficient: 0.000000000001))
	
	public static let tesseractFeet = UnitTesseract(symbol: "ft⁴", converter: UnitConverterLinear(coefficient: 0.00863097485))
	public static let tesseractInches = UnitTesseract(symbol: "in⁴", converter: UnitConverterLinear(coefficient: 0.000000416231426))
	
	public override class func baseUnit() -> Self {
		return UnitTesseract.tesseractMeters as! Self
	}
	public static let allEngineeringUnits: [UnitTesseract] = [.tesseractMeters,.tesseractCentimeters,.tesseractMillimeters,.tesseractFeet,.tesseractInches]
	override public static var supportsSecureCoding: Bool { true }
	public var isImperial: Bool {
		if self ==|| [.tesseractFeet,.tesseractInches] {
			return true
		} else { return false }
	}
	/// The unit that results by taking the square root of a UnitTesseract unit.
	public var squareBaseUnit: UnitArea! {
		switch self {
		case .tesseractMeters: 		return UnitArea.squareMeters
		case .tesseractCentimeters: return UnitArea.squareCentimeters
		case .tesseractMillimeters: return UnitArea.squareMillimeters
		case .tesseractFeet: 		return UnitArea.squareFeet
		case .tesseractInches: 		return UnitArea.squareInches
		default:					return nil
		}
	}
}
public extension Measurement where UnitType: UnitTesseract {
	func sqrt() -> Measurement<UnitArea> {
		if let squareUnit = self.unit.squareBaseUnit {
			return Measurement<UnitArea>(value: Darwin.sqrt(self.value), unit: squareUnit)
		} else {
			return Measurement<UnitArea>(value: Darwin.sqrt(self.baseUnitValue()), unit: .baseUnit())
		}
	}
}
extension UnitAcceleration: EngineeringUnit {
	public static let feetPerSecondSquared = UnitAcceleration(symbol: "ft/s²", converter: UnitConverterLinear(coefficient: 3.2808398950131))
	public static let imperialGravity = UnitAcceleration(symbol: "g (ft/s²)", converter: UnitConverterLinear(coefficient: 32.18503937007874))
	public static let allEngineeringUnits: [UnitAcceleration] = [.metersPerSecondSquared,.gravity,.feetPerSecondSquared,.imperialGravity]
	public var isImperial: Bool {
		if self ==|| [.feetPerSecondSquared,.imperialGravity] {
			return true
		} else { return false }
	}
}
extension UnitPressure: EngineeringUnit {
	public static let kipsPerSquareInch = UnitPressure(symbol: "ksi", converter: UnitConverterLinear(coefficient: 6894757.2932))
	public static let poundsForcePerSquareFoot = UnitPressure(symbol: "psf", converter: UnitConverterLinear(coefficient: 47.88025898))
	public static let allEngineeringUnits: [UnitPressure] = [.gigapascals,.megapascals,.kilopascals,.newtonsPerMetersSquared,.kipsPerSquareInch,.poundsForcePerSquareFoot,.poundsForcePerSquareInch]
	public var isImperial: Bool {
		if self ==|| [.inchesOfMercury,.poundsForcePerSquareInch,.kipsPerSquareInch,.poundsForcePerSquareFoot,.bars,.millibars] {
			return true
		} else { return false }
	}
}
extension UnitTemperature: EngineeringUnit {
	public static let allEngineeringUnits: [UnitTemperature] = [.celsius,.fahrenheit,.kelvin]
	
	public var positiveOnly: Bool { self == .kelvin }
	
	public var isImperial: Bool {
		if self == .fahrenheit {
			return true
		} else { return false }
	}
}
extension UnitAngle: EngineeringUnit {
	public var isImperial: Bool { self == UnitAngle.degrees }
	public static let allEngineeringUnits: [UnitAngle] = [.radians,.gradians,.degrees,.revolutions]
	public static var allImperialEngineeringUnitSymbols: [String] { allEngineeringUnits.map({$0.symbol}) }
	public static var allSIEngineeringUnitSymbols: [String] { allEngineeringUnits.map({$0.symbol}) }
}

public extension Measurement {
	
	// ρ = m/V
	static func * <R: UnitMass>(a: Measurement<UnitVolume>, b: Measurement<UnitDensity>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitMass>(a: Measurement<UnitDensity>, b: Measurement<UnitVolume>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitDensity>(a: Measurement<UnitMass>, b: Measurement<UnitVolume>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitVolume>(a: Measurement<UnitMass>, b: Measurement<UnitDensity>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// M = Fd
	static func * <R: UnitWork>(a: Measurement<UnitForce>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitWork>(a: Measurement<UnitLength>, b: Measurement<UnitForce>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitForce>(a: Measurement<UnitWork>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitWork>, b: Measurement<UnitForce>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// F = W*d
	static func * <R: UnitForce>(a: Measurement<UnitLinearForce>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement<UnitLength>, b: Measurement<UnitLinearForce>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLinearForce>(a: Measurement<UnitForce>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitForce>, b: Measurement<UnitLinearForce>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// F/L = P*L
	static func * <R: UnitLinearForce>(a: Measurement<UnitPressure>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitLinearForce>(a: Measurement<UnitLength>, b: Measurement<UnitPressure>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitPressure>(a: Measurement<UnitLinearForce>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitLinearForce>, b: Measurement<UnitPressure>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// F = ma
	static func * <R: UnitForce>(a: Measurement<UnitMass>, b: Measurement<UnitAcceleration>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement<UnitAcceleration>, b: Measurement<UnitMass>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitAcceleration>(a: Measurement<UnitForce>, b: Measurement<UnitMass>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitMass>(a: Measurement<UnitForce>, b: Measurement<UnitAcceleration>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// P = F/A
	static func * <R: UnitForce>(a: Measurement<UnitArea>, b: Measurement<UnitPressure>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement<UnitPressure>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitPressure>(a: Measurement<UnitForce>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitArea>(a: Measurement<UnitForce>, b: Measurement<UnitPressure>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// Tesseract
	static func * <R: UnitTesseract>(a: Measurement<UnitVolume>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitTesseract>(a: Measurement<UnitLength>, b: Measurement<UnitVolume>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitTesseract>(a: Measurement<UnitArea>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitTesseract>, b: Measurement<UnitVolume>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitVolume>(a: Measurement<UnitTesseract>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitArea>(a: Measurement<UnitTesseract>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// V = A*L
	static func * <R: UnitVolume>(a: Measurement<UnitArea>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitVolume>(a: Measurement<UnitLength>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitVolume>, b: Measurement<UnitArea>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitArea>(a: Measurement<UnitVolume>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	// A = L*L
	static func * <R: UnitArea>(a: Measurement<UnitLength>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement<UnitArea>, b: Measurement<UnitLength>) -> Measurement<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// T*(1/T) = 1
	// Note that Double/Measurement always returns with the same unit as the Measurement!
	static func * (a: Measurement<UnitTemperature>, b: Measurement<UnitInverseTemperature>) -> Double {
		return a.converted(to: .baseUnit()).value * b.converted(to: .baseUnit()).value
	}
	static func * (a: Measurement<UnitInverseTemperature>, b: Measurement<UnitTemperature>) -> Double {
		return a.converted(to: .baseUnit()).value * b.converted(to: .baseUnit()).value
	}
	
	/// A function that produces a dimensionless scalar from two units with the same `Dimension` type.
	static func / <D: Dimension>(lhs: Measurement<D>, rhs: Measurement<D>) -> Double {
		return lhs.converted(to: .baseUnit()).value / rhs.converted(to: .baseUnit()).value
	}
	
	// Private function to create multiplicative relationships between measurements
	private static func multiplyMeasurementsIntoNewMeasurement<A: Dimension, B: Dimension, R: Dimension>(a: Measurement<A>, b: Measurement<B>, returnType: R.Type) -> Measurement<R> {
		let value = a.converted(to: .baseUnit()).value * b.converted(to: .baseUnit()).value
		return Measurement<R>(value: value, unit: .baseUnit())
	}
	// Private function to create divisor relationships between measurements
	private static func divideMeasurementsIntoNewMeasurement<A: Dimension, B: Dimension, R: Dimension>(a: Measurement<A>, b: Measurement<B>, returnType: R.Type) -> Measurement<R> {
		let value = a.converted(to: .baseUnit()).value / b.converted(to: .baseUnit()).value
		return Measurement<R>(value: value, unit: .baseUnit())
	}
	// Negation
	static prefix func - <D: EngineeringUnit>(measurement: Measurement<D>) -> Measurement<D> {
		if measurement.unit.positiveOnly {
			return measurement
		} else {
			return Measurement<D>(value: -measurement.value, unit: measurement.unit)
		}
	}
}

public extension Measurement where UnitType: Dimension {
	init() {
		self.init(value: 0, unit: .baseUnit())
	}
	func baseUnitValue() -> Double {
		return self.converted(to: .baseUnit()).value
	}
	/// Returns the absolute value of a Measurement.
	func abs() -> Self {
		Measurement(value: Swift.abs(self.value), unit: self.unit)
	}
}
extension Measurement: @retroactive AdditiveArithmetic where UnitType: Dimension {
	public static var zero: Measurement<UnitType> {
		Self()
	}
}
//extension Measurement.FormatStyle: ParseableFormatStyle {
//	public typealias Strategy = MeasurementFormatStyleParsingStrategy
//	public var parseStrategy: Self.Strategy<UnitType> {
//		MeasurementFormatStyleParsingStrategy()
//	}
//}
//// where Self.FormatInput == Self.Strategy.ParseOutput, Self.FormatOutput == Self.Strategy.ParseInput
//public struct MeasurementFormatStyleParsingStrategy<UnitType: Dimension>: ParseStrategy {
//	public func parse(_ value: String) throws -> Measurement<UnitType> {
//		value.formatted(.measurement(width: .abbreviated, usage: .asProvided, numberFormatStyle: .localizedDouble(locale: Locale.current)))
//	}
//}


import simd
import SwiftUI
public struct Measurement3D<UnitType>: Comparable, CustomDebugStringConvertible, CustomStringConvertible, Hashable, Codable where UnitType: Unit {
	
	/// The unit component of the measurement.
	@CodableViaNSCoding public var unit: UnitType
	/// The value components of the measurement.
	public var values: SIMD3<Double>
	/// The x-value component of the measurement.
	public var xMeasurement: Measurement<UnitType> { Measurement(value: values.x, unit: unit) }
	/// The y-value component of the measurement.
	public var yMeasurement: Measurement<UnitType> { Measurement(value: values.y, unit: unit) }
	/// The z-value component of the measurement.
	public var zMeasurement: Measurement<UnitType> { Measurement(value: values.z, unit: unit) }
	
	public init(values: SIMD3<Double> = SIMD3(), unit: UnitType) {
		_unit = CodableViaNSCoding(wrappedValue: unit)
		self.values = values
	}
	
	
	
	public var description: String {
		"(x: \(xMeasurement), y: \(yMeasurement), z: \(zMeasurement))"
	}
	public var debugDescription: String {
		"""
		unit: \(unit.debugDescription)
		x: \(xMeasurement.debugDescription)
		y: \(yMeasurement.debugDescription)
		z: \(zMeasurement.debugDescription)
		"""
	}
	
	public static func < (lhs: Measurement3D<UnitType>, rhs: Measurement3D<UnitType>) -> Bool {
		simd_length(lhs.values) < simd_length(rhs.values)
	}
	
	/// Multiply a scalar value by a measurement.
	public static func * (lhs: Double, rhs: Measurement3D<UnitType>) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs*rhs.values, unit: rhs.unit)
	}
	/// Multiply a measurement by a scalar value.
	public static func * (lhs: Measurement3D<UnitType>, rhs: Double) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs.values*rhs, unit: lhs.unit)
	}
	/// Adds one measurement to another.
	public static func + (lhs: Measurement3D<UnitType>, rhs: Measurement3D<UnitType>) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs.values+rhs.values, unit: lhs.unit)
	}
	/// Subtract one measurement from another.
	public static func - (lhs: Measurement3D<UnitType>, rhs: Measurement3D<UnitType>) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs.values-rhs.values, unit: lhs.unit)
	}
	/// Divide a scalar value by a measurement.
	public static func / (lhs: Double, rhs: Measurement3D<UnitType>) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs/rhs.values, unit: rhs.unit)
	}
	/// Divide a measurement by a scalar value.
	public static func / (lhs: Measurement3D<UnitType>, rhs: Double) -> Measurement3D<UnitType> {
		Measurement3D<UnitType>(values: lhs.values/rhs, unit: lhs.unit)
	}
}
public extension Measurement3D where UnitType: Dimension {
	
	func baseUnitValues() -> SIMD3<Double> {
		return self.converted(to: .baseUnit()).values
	}
	
	/// Converts the measurement to the specified unit.
	mutating func convert(to newUnit: UnitType) {
		values.x = xMeasurement.converted(to: newUnit).value
		values.y = yMeasurement.converted(to: newUnit).value
		values.z = zMeasurement.converted(to: newUnit).value
		unit = newUnit
	}
	/// Returns a new measurement created by converting to the specified unit.
	func converted(to newUnit: UnitType) -> Measurement3D<UnitType> {
		var copy = self
		copy.convert(to: newUnit)
		return copy
	}
	/// A function that produces a dimensionless vector from two units with the same `Dimension` type.
	static func / <D: Dimension>(lhs: Measurement3D<D>, rhs: Measurement3D<D>) -> SIMD3<Double> {
		let v0 = lhs.converted(to: .baseUnit()).values.x / rhs.converted(to: .baseUnit()).values.x
		let v1 = lhs.converted(to: .baseUnit()).values.y / rhs.converted(to: .baseUnit()).values.y
		let v2 = lhs.converted(to: .baseUnit()).values.z / rhs.converted(to: .baseUnit()).values.z
		return SIMD3(v0, v1, v2)
	}
	/// Returns the effective "distance" of the measurement using pythagorean theorem.
	func distance() -> Measurement<UnitType> {
		Measurement(value: sqrt(values.x**2+values.y**2+values.z**2), unit: self.unit)
	}
	init() {
		self.init(values: SIMD3(), unit: .baseUnit())
	}
}
extension Measurement3D: AdditiveArithmetic where UnitType: Dimension {
	public static var zero: Measurement3D<UnitType> {
		Self()
	}
}
public extension Measurement3D {
	
	// M = Fd
	static func * <R: UnitWork>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitWork>(a: Measurement3D<UnitLength>, b: Measurement3D<UnitForce>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitForce>(a: Measurement3D<UnitWork>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement3D<UnitWork>, b: Measurement3D<UnitForce>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitWork>(a: Measurement3D<UnitForce>, b: Measurement<UnitLength>) -> Measurement3D<R> {
		let b2 = Measurement3D<UnitLength>(values: SIMD3(b.value,b.value,b.value), unit: b.unit)
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b2, returnType: R.self)
	}
	static func * <R: UnitWork>(a: Measurement<UnitLength>, b: Measurement3D<UnitForce>) -> Measurement3D<R> {
		let a2 = Measurement3D<UnitLength>(values: SIMD3(a.value,a.value,a.value), unit: a.unit)
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a2, b: b, returnType: R.self)
	}
	static func / <R: UnitForce>(a: Measurement3D<UnitWork>, b: Measurement<UnitLength>) -> Measurement3D<R> {
		let b2 = Measurement3D<UnitLength>(values: SIMD3(b.value,b.value,b.value), unit: b.unit)
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b2, returnType: R.self)
	}
	
	// F = W*d
	static func * <R: UnitForce>(a: Measurement3D<UnitLinearForce>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement3D<UnitLength>, b: Measurement3D<UnitLinearForce>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLinearForce>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitLinearForce>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement3D<UnitLinearForce>, b: Measurement<UnitLength>) -> Measurement3D<R> {
		let b2 = Measurement3D<UnitLength>(values: SIMD3(b.value,b.value,b.value), unit: b.unit)
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b2, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement<UnitLength>, b: Measurement3D<UnitLinearForce>) -> Measurement3D<R> {
		let a2 = Measurement3D<UnitLength>(values: SIMD3(a.value,a.value,a.value), unit: a.unit)
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a2, b: b, returnType: R.self)
	}
	static func / <R: UnitLinearForce>(a: Measurement3D<UnitForce>, b: Measurement<UnitLength>) -> Measurement3D<R> {
		let b2 = Measurement3D<UnitLength>(values: SIMD3(b.value,b.value,b.value), unit: b.unit)
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b2, returnType: R.self)
	}
	
	// W = P*L
	static func * <R: UnitLinearForce>(a: Measurement3D<UnitPressure>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitLinearForce>(a: Measurement3D<UnitLength>, b: Measurement3D<UnitPressure>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitPressure>(a: Measurement3D<UnitLinearForce>, b: Measurement3D<UnitLength>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitLength>(a: Measurement3D<UnitLinearForce>, b: Measurement3D<UnitPressure>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// F = ma
	static func * <R: UnitForce>(a: Measurement3D<UnitMass>, b: Measurement3D<UnitAcceleration>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement3D<UnitAcceleration>, b: Measurement3D<UnitMass>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitAcceleration>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitMass>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitMass>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitAcceleration>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	// P = F/A
	static func * <R: UnitForce>(a: Measurement3D<UnitArea>, b: Measurement3D<UnitPressure>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func * <R: UnitForce>(a: Measurement3D<UnitPressure>, b: Measurement3D<UnitArea>) -> Measurement3D<R> {
		return Self.multiplyMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitPressure>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitArea>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	static func / <R: UnitArea>(a: Measurement3D<UnitForce>, b: Measurement3D<UnitPressure>) -> Measurement3D<R> {
		return Self.divideMeasurementsIntoNewMeasurement(a: a, b: b, returnType: R.self)
	}
	
	/// A function that produces a dimensionless scalar from two units with the same `Dimension` type.
	static func / <D: Dimension>(lhs: Measurement3D<D>, rhs: Measurement3D<D>) -> SIMD3<Double> {
		return lhs.converted(to: .baseUnit()).values / rhs.converted(to: .baseUnit()).values
	}
	
	// Private function to create multiplicative relationships between measurements
	private static func multiplyMeasurementsIntoNewMeasurement<A: Dimension, B: Dimension, R: Dimension>(a: Measurement3D<A>, b: Measurement3D<B>, returnType: R.Type) -> Measurement3D<R> {
		let values = a.converted(to: .baseUnit()).values * b.converted(to: .baseUnit()).values
		return Measurement3D<R>(values: values, unit: .baseUnit())
	}
	// Private function to create divisor relationships between measurements
	private static func divideMeasurementsIntoNewMeasurement<A: Dimension, B: Dimension, R: Dimension>(a: Measurement3D<A>, b: Measurement3D<B>, returnType: R.Type) -> Measurement3D<R> {
		let values = a.converted(to: .baseUnit()).values / b.converted(to: .baseUnit()).values
		return Measurement3D<R>(values: values, unit: .baseUnit())
	}
}

public enum UnitSystem: String, CaseIterable, Identifiable, Hashable {
	case imperial = "Imperial", SI = "SI"
	public var id: String { rawValue }
	public static func selection(for storedData: String) -> [UnitSystem] {
		if storedData == "all" { return UnitSystem.allCases }
		if storedData == "Imperial" { return [.imperial] }
		if storedData == "SI" { return [.SI] }
		else { return [.imperial] }
	}
}
public extension Array where Element == UnitSystem {
	var storedDataString: String {
		if self.count > 1 { "all" } else if self.count == 1 { self.first!.rawValue } else { "Imperial" }
	}
}

public struct ENGRValueField<EngrUnitType: EngineeringUnit>: View where EngrUnitType == EngrUnitType.EngDimension {
	
	@Environment(\.deviceOS) var os
	let description: String
	@Binding var measurement: Measurement<EngrUnitType>
	@State private var measurementUnit: String
	var measurementValue: Binding<Double> {
		Binding {
			measurement.converted(to: getUnit()).value
		} set: { newValue in
			measurement = Measurement(value: newValue, unit: getUnit())
		}
	}
	//let minValue: Measurement<EngrUnitType>?
	//let maxValue: Measurement<EngrUnitType>?
	@FocusState var thisMeasurementIsFocused: Bool
	let positiveOnly: Bool
	@AppStorage("preferredUnitSystem") var preferredUnitsData: String = "Imperial"
	@State private var allowedUnitSystems: [UnitSystem]
	let fixedUnitSystem: Bool
	let defaultImperialUnitSymbol: String
	let defaultSIUnitSymbol: String
	
	public init(_ description: String, _ measurement: Binding<Measurement<EngrUnitType>>, allowedUnits: [UnitSystem], positiveOnly: Bool = false)  {
		self.description = description
		_measurement = measurement
		_measurementUnit = State(initialValue: measurement.wrappedValue.unit.symbol)
		self.positiveOnly = positiveOnly
		_allowedUnitSystems = State(initialValue: allowedUnits)
		self.fixedUnitSystem = true
		defaultImperialUnitSymbol = EngrUnitType.allImperialEngineeringUnitSymbols.first!
		defaultSIUnitSymbol = EngrUnitType.allSIEngineeringUnitSymbols.first!
	}
//	public init(_ description: String, _ measurement: Binding<Measurement<EngrUnitType>>, allowedUnits: [UnitSystem]? = nil, minValue: Measurement<EngrUnitType>? = nil, maxValue: Measurement<EngrUnitType>? = nil, positiveOnly: Bool = false)  {
//		self.description = description
//		_measurement = measurement
//		self.minValue = minValue
//		self.maxValue = maxValue
//		_measurementUnit = State(initialValue: measurement.wrappedValue.unit.symbol)
//		self.positiveOnly = positiveOnly
//		self.allowedUnitSystems = allowedUnits ?? UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial")
//	}
	public init(_ description: String, _ measurement: Binding<Measurement<EngrUnitType>>, defaultImperialUnit: EngrUnitType, defaultSIUnit: EngrUnitType, positiveOnly: Bool = false)  {
		self.description = description
		_measurement = measurement
		_measurementUnit = State(initialValue: measurement.wrappedValue.unit.symbol)
		self.positiveOnly = positiveOnly
		_allowedUnitSystems = State(initialValue: UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial"))
		self.fixedUnitSystem = false
		self.defaultImperialUnitSymbol = defaultImperialUnit.symbol
		self.defaultSIUnitSymbol = defaultSIUnit.symbol
	}
	
	let measurementFormatStyle: Measurement<EngrUnitType>.FormatStyle = .measurement(width: .abbreviated, usage: .asProvided, numberFormatStyle: .localizedDouble(locale: Locale.current))
	
	public var body: some View {
		HStack {
			Text("\(description)")
			Spacer()
			TextField(description, value: measurementValue, format: FloatingPointMathParseableFormatStyle(), prompt: Text(""))
				.textFieldStyle(.roundedBorder)
#if !os(macOS)
				.keyboardType(.numbersAndPunctuation)
				.frame(minWidth: 80, idealWidth: 100, maxWidth: 140)
#else
				.frame(minWidth: 80, idealWidth: 100, maxWidth: 120)
#endif
				.focused($thisMeasurementIsFocused)
				.onChange(of: measurement) {
					if measurement.value < 0 && positiveOnly {
						measurement = Measurement(value: abs(measurement.value), unit: measurement.unit)
					}
					measurementUnit = measurement.unit.symbol
				}
			Picker(os == .macOS ? "":"Unit for \(description)", selection: $measurementUnit) {
				if allowedUnitSystems.contains(.imperial) {
					ForEach(EngrUnitType.allImperialEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
				if allowedUnitSystems.contains(.SI) {
					ForEach(EngrUnitType.allSIEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
			}
			.onChange(of: measurementUnit) {
				let unit = getUnit()
				measurement = Measurement(value: measurement.value, unit: unit) // Just change units
				//measurement = Measurement(value: measurement.converted(to: unit).value, unit: unit) // Conversion Option
			}
			.macOS { $0.frame(minWidth: 70, idealWidth: 100, maxWidth: 120) }
		}
		.onChange(of: preferredUnitsData) { oldValue, newValue in
			if !fixedUnitSystem {
				self.allowedUnitSystems = UnitSystem.selection(for: preferredUnitsData)
				if allowedUnitSystems == [.imperial] {
					measurementUnit = defaultImperialUnitSymbol
				}
				if allowedUnitSystems == [.SI] {
					measurementUnit = defaultSIUnitSymbol
				}
			}
		}
//		.onAppear {
//			if !fixedUnitSystem {
//				self.allowedUnitSystems = UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial")
//				let isImperial = EngrUnitType.allImperialEngineeringUnitSymbols.contains(measurementUnit)
//				if isImperial && !allowedUnitSystems.contains(.imperial) || !isImperial && !allowedUnitSystems.contains(.SI) {
//					if allowedUnitSystems == [.imperial] {
//						measurementUnit = defaultImperialUnitSymbol
//					}
//					if allowedUnitSystems == [.SI] {
//						measurementUnit = defaultSIUnitSymbol
//					}
//				}
//			}
//		}
//		.toolbar {
//			if thisMeasurementIsFocused {
//				ToolbarItemGroup(placement: .keyboard) {
//					if !measurement.unit.positiveOnly && !positiveOnly {
//						Button("Negate", systemImage: "plus.forwardslash.minus") {
//							measurement = Measurement(value: -measurement.value, unit: measurement.unit)
//						}
//					}
//					Spacer()
//					Button {
//						thisMeasurementIsFocused = false
//					} label: {
//						Text("Done").foregroundStyle(Color.accentColor)
//					}
//				}
//			}
//		}
	}
	
	func getUnit() -> EngrUnitType {
		let unitIndex = EngrUnitType.allEngineeringUnitSymbols.firstIndex(of: measurementUnit)!
		return EngrUnitType.allEngineeringUnits[unitIndex]
	}
}
public struct ENGRValueDisplay<EngrUnitType: EngineeringUnit>: View where EngrUnitType == EngrUnitType.EngDimension {
	
	@Environment(\.deviceOS) var os
	let description: String
	@State var measurement: Measurement<EngrUnitType>
	@State private var measurementUnit: String
	@AppStorage("preferredUnitSystem") var preferredUnitsData: String = "Imperial"
	@State private var allowedUnitSystems: [UnitSystem]
	let fixedUnitSystem: Bool
	let defaultImperialUnitSymbol: String
	let defaultSIUnitSymbol: String
	
	public init(_ description: String, _ measurement: Measurement<EngrUnitType>, allowedUnits: [UnitSystem])  {
		self.description = description
		_measurement = State(initialValue: measurement)
		_measurementUnit = State(initialValue: measurement.unit.symbol)
		_allowedUnitSystems = State(initialValue: allowedUnits)
		self.fixedUnitSystem = true
		defaultImperialUnitSymbol = EngrUnitType.allImperialEngineeringUnitSymbols.first!
		defaultSIUnitSymbol = EngrUnitType.allSIEngineeringUnitSymbols.first!
	}
	public init(_ description: String, _ measurement: Measurement<EngrUnitType>, defaultImperialUnit: EngrUnitType, defaultSIUnit: EngrUnitType)  {
		self.description = description
		_measurement = State(initialValue: measurement)
		_measurementUnit = State(initialValue: measurement.unit.symbol)
		_allowedUnitSystems = State(initialValue: UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial"))
		self.fixedUnitSystem = false
		self.defaultImperialUnitSymbol = defaultImperialUnit.symbol
		self.defaultSIUnitSymbol = defaultSIUnit.symbol
	}
	
	let measurementFormatStyle: Measurement<EngrUnitType>.FormatStyle = .measurement(width: .abbreviated, usage: .asProvided, numberFormatStyle: .localizedDouble(locale: Locale.current))
	
	public var body: some View {
		HStack {
			Text(description)
			Spacer()
			Text(measurement.converted(to: getUnit()).value.formatted(sigFigs: ...4))
			Picker(os == .macOS ? "":"Unit for \(description)", selection: $measurementUnit) {
				if allowedUnitSystems.contains(.imperial) {
					ForEach(EngrUnitType.allImperialEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
				if allowedUnitSystems.contains(.SI) {
					ForEach(EngrUnitType.allSIEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
			}
			.onChange(of: measurementUnit) {
				let unit = getUnit()
				measurement = Measurement(value: measurement.converted(to: unit).value, unit: unit)
			}
			.pickerStyle(.menu)
			.macOS { $0.frame(minWidth: 70, idealWidth: 100, maxWidth: 120) }
		}
		.onChange(of: preferredUnitsData) { oldValue, newValue in
			if !fixedUnitSystem {
				self.allowedUnitSystems = UnitSystem.selection(for: preferredUnitsData)
				if allowedUnitSystems == [.imperial] {
					measurementUnit = defaultImperialUnitSymbol
				}
				if allowedUnitSystems == [.SI] {
					measurementUnit = defaultSIUnitSymbol
				}
			}
		}
//		.onAppear {
//			if !fixedUnitSystem {
//				self.allowedUnitSystems = UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial")
//				let isImperial = EngrUnitType.allImperialEngineeringUnitSymbols.contains(measurementUnit)
//				if isImperial && !allowedUnitSystems.contains(.imperial) || !isImperial && !allowedUnitSystems.contains(.SI) {
//					if allowedUnitSystems == [.imperial] {
//						measurementUnit = defaultImperialUnitSymbol
//					}
//					if allowedUnitSystems == [.SI] {
//						measurementUnit = defaultSIUnitSymbol
//					}
//				}
//			}
//		}
	}
	
	func getUnit() -> EngrUnitType {
		let unitIndex = EngrUnitType.allEngineeringUnitSymbols.firstIndex(of: measurementUnit)!
		return EngrUnitType.allEngineeringUnits[unitIndex]
	}
}
public struct ENGRMeasurementPicker<EngrUnitType: EngineeringUnit>: View where EngrUnitType == EngrUnitType.EngDimension {
	
	let description: String
	@Binding var unit: EngrUnitType
	@State private var measurementUnit: String
	@AppStorage("preferredUnitSystem") var preferredUnitsData: String = "Imperial"
	@State private var allowedUnitSystems: [UnitSystem]
	let fixedUnitSystem: Bool
	let defaultImperialUnitSymbol: String
	let defaultSIUnitSymbol: String
	
	public init(_ description: String, unit: Binding<EngrUnitType>, allowedUnits: [UnitSystem])  {
#if os(macOS)
		self.description = description+":"
#else
		self.description = description
#endif
		_unit = unit
		_measurementUnit = State(initialValue: unit.wrappedValue.symbol)
		_allowedUnitSystems = State(initialValue: allowedUnits)
		self.fixedUnitSystem = true
		defaultImperialUnitSymbol = EngrUnitType.allImperialEngineeringUnitSymbols.first!
		defaultSIUnitSymbol = EngrUnitType.allSIEngineeringUnitSymbols.first!
	}
	public init(_ description: String, unit: Binding<EngrUnitType>, defaultImperialUnit: EngrUnitType, defaultSIUnit: EngrUnitType)  {
#if os(macOS)
		self.description = description+":"
#else
		self.description = description
#endif
		_unit = unit
		_measurementUnit = State(initialValue: unit.wrappedValue.symbol)
		_allowedUnitSystems = State(initialValue: UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial"))
		self.fixedUnitSystem = false
		self.defaultImperialUnitSymbol = defaultImperialUnit.symbol
		self.defaultSIUnitSymbol = defaultSIUnit.symbol
	}
	
	let measurementFormatStyle: Measurement<EngrUnitType>.FormatStyle = .measurement(width: .abbreviated, usage: .asProvided, numberFormatStyle: .localizedDouble(locale: Locale.current))
	
	public var body: some View {
		HStack {
//			Text(description)
//			Spacer()
			Picker(description, selection: $measurementUnit) {
				if allowedUnitSystems.contains(.imperial) {
					ForEach(EngrUnitType.allImperialEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
				if allowedUnitSystems.contains(.SI) {
					ForEach(EngrUnitType.allSIEngineeringUnitSymbols, id: \.self) { unitSymbol in
						Text(unitSymbol).tag(unitSymbol)
					}
				}
			}
			.onChange(of: measurementUnit) {
				unit = getUnit()
			}
			.pickerStyle(.menu)
		}
		.onChange(of: preferredUnitsData) { oldValue, newValue in
			if !fixedUnitSystem {
				self.allowedUnitSystems = UnitSystem.selection(for: preferredUnitsData)
				if allowedUnitSystems == [.imperial] {
					measurementUnit = defaultImperialUnitSymbol
				}
				if allowedUnitSystems == [.SI] {
					measurementUnit = defaultSIUnitSymbol
				}
			}
		}
//		.onAppear {
//			if !fixedUnitSystem {
//				self.allowedUnitSystems = UnitSystem.selection(for: UserDefaults.standard.string(forKey: "preferredUnitSystem") ?? "Imperial")
//				let isImperial = EngrUnitType.allImperialEngineeringUnitSymbols.contains(measurementUnit)
//				if isImperial && !allowedUnitSystems.contains(.imperial) || !isImperial && !allowedUnitSystems.contains(.SI) {
//					if allowedUnitSystems == [.imperial] {
//						measurementUnit = defaultImperialUnitSymbol
//					}
//					if allowedUnitSystems == [.SI] {
//						measurementUnit = defaultSIUnitSymbol
//					}
//				}
//			}
//		}
	}
	
	func getUnit() -> EngrUnitType {
		let unitIndex = EngrUnitType.allEngineeringUnitSymbols.firstIndex(of: measurementUnit)!
		return EngrUnitType.allEngineeringUnits[unitIndex]
	}
}

struct FieldsPreviews: PreviewProvider {
	static var previews: some View {
		ZStack {
			VStack {
				ENGRValueField("Length", .constant(Measurement<UnitDensity>(value: 12, unit: .kilogramPerCubicMeter)), allowedUnits: [.imperial,.SI])
					.padding()
				Spacer()
				ENGRValueDisplay("Length", Measurement<UnitDensity>(value: 12.0110010, unit: .kilogramPerCubicMeter), allowedUnits: [.imperial,.SI])
					.padding()
			}
		}
	}
}

public struct FloatingPointMathParseableFormatStyle: ParseableFormatStyle {
	
	public var parseStrategy = MathParser()
	public func format(_ value: Double) -> String {
		value.formatted(.number)
	}
	
	public struct MathParser: ParseStrategy {
		public func parse(_ value: String) throws -> Double {
//			let expn = NSExpression(format:value)
//			if let result = expn.expressionValue(with: nil, context: nil) as? Double {
//				return result
//			} else {
//				print("Can't evaluate \(value) with NSExpression.")
//				return try Double(value, format: .number)
//			}
			do {
				let expression = try MathExpression(value)
				return expression.evaluate()
			} catch {
				return try Double(value, format: .number)
			}
		}
	}
}



