import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_COLUMNS = ['X [m]', 'Y [m]', 'Z[m]', 'Ex [V/m]', 'Ey [V/m]', 'Ez [V/m]']


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="Plot Ey vs X for a given Y slice from ANSYS field data.")
	parser.add_argument(
		"input_file",
		nargs="?",
		help="Path to the ANSYS field export (whitespace-delimited). If omitted, the script looks for 'ansys_field.txt'.",
	)
	parser.add_argument(
		"--y-target",
		type=float,
		default=12e-6,
		help="Target Y position in meters (default: 12e-6).",
	)
	parser.add_argument(
		"--skiprows",
		type=int,
		default=0,
		help="Number of leading rows to skip when reading the file (default: 0).",
	)
	parser.add_argument(
		"--delimiter",
		default=None,
		help="Custom delimiter to use instead of whitespace (default: whitespace).",
	)
	parser.add_argument(
		"--output",
		default="Ey_vs_X.png",
		help="Path of the PNG file to save (default: Ey_vs_X.png).",
	)
	parser.add_argument(
		"--show",
		action="store_true",
		help="Display the plot window in addition to saving the PNG.",
	)
	return parser.parse_args()


def resolve_input_path(arg_value: str | None) -> Path:
	if arg_value:
		candidate = Path(arg_value).expanduser()
		if candidate.is_dir():
			raise ValueError(f"'{candidate}' è una cartella, serve un file dati.")
		return candidate

	# Fallback to a sensible default next to the script
	default_path = Path(__file__).with_suffix('.txt')
	if default_path.exists():
		return default_path

	raise FileNotFoundError(
		"Nessun file di input specificato e non esiste 'ansys_field.txt'. "
		"Passa il percorso del file dati come argomento."
	)


def load_field_dataframe(
	file_path: Path,
	skiprows: int,
	delimiter: str | None,
) -> pd.DataFrame:
	if file_path.suffix == '.py':
		raise ValueError(
			"Il file di input punta allo script Python stesso. Specifica il file dati esportato da ANSYS."
		)

	read_kwargs: dict[str, object] = {
		"delim_whitespace": delimiter is None,
		"on_bad_lines": "skip",
	}
	if delimiter is not None:
		read_kwargs["sep"] = delimiter

	df = pd.read_csv(
		file_path,
		names=DEFAULT_COLUMNS,
		skiprows=skiprows,
		header=None,
		**read_kwargs,
	)

	# Coerce columns to numeric values; non-numeric rows become NaN and can be dropped.
	for col in DEFAULT_COLUMNS:
		if col in df.columns:
			df[col] = pd.to_numeric(df[col], errors='coerce')

	# Drop rows lacking the essential numeric information
	df = df.dropna(subset=['X [m]', 'Y [m]', 'Ey [V/m]']).reset_index(drop=True)

	if df.empty:
		raise ValueError(
			"Nessun dato numerico trovato dopo la pulizia. Controlla 'skiprows' e il formato del file."
		)

	return df


def find_nearest_y(df: pd.DataFrame, target_y: float) -> float:
	y_values = df['Y [m]'].to_numpy()
	nearest_idx = np.abs(y_values - target_y).argmin()
	return float(y_values[nearest_idx])


def main() -> None:
	args = parse_args()
	input_path = resolve_input_path(args.input_file)

	df = load_field_dataframe(input_path, skiprows=args.skiprows, delimiter=args.delimiter)

	y_nearest = find_nearest_y(df, args.y_target)
	print(f"Y target: {args.y_target:.3e} m, Y più vicino nel file: {y_nearest:.3e} m")

	df_slice = df[np.isclose(df['Y [m]'], y_nearest, atol=1e-12)]
	if df_slice.empty:
		# As a fallback, use tolerance-based selection
		tolerance = np.abs(y_nearest - args.y_target) * 1.01 + 1e-12
		df_slice = df[np.abs(df['Y [m]'] - y_nearest) <= tolerance]

	if df_slice.empty:
		raise ValueError("Impossibile trovare dati per il valore di Y selezionato.")

	plt.figure(figsize=(10, 6))
	plt.plot(df_slice['X [m]'], df_slice['Ey [V/m]'], label=f'Y = {y_nearest*1e6:.2f} µm')
	plt.title('Campo Elettrico Ey lungo X a Y specifico')
	plt.xlabel('X [m]')
	plt.ylabel('Ey [V/m]')
	plt.grid(True)
	plt.legend()
	plt.tight_layout()
	plt.savefig(args.output, dpi=300)
	print(f"Grafico salvato in '{args.output}'.")

	if args.show:
		plt.show()


if __name__ == '__main__':
	main()

