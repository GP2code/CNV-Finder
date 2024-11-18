import argparse
import numpy as np
from cnv_finder.model_methods import prep_ml_datasets, train_binary_lstm, model_predict


def main():
    parser = argparse.ArgumentParser(
        description='Arguments for Running CNV-Finder ML Model Training/Prediction.')
    parser.add_argument('--train_file', type=str, default=None,
                        help='Path to training file with windowed samples.')
    parser.add_argument('--test_file', type=str, default=None,
                        help='Path to testing file with windowed samples.')
    parser.add_argument('--train', action='store_true',
                        help='Train new model when True.')
    parser.add_argument('--feature_names', type=str, nargs='*',
                        help='Add feature names for model training.')
    parser.add_argument('--model_file', type=str, default=None,
                        help='Path to saved model file (.keras).')  # stil need to fix keras pkl issue
    parser.add_argument('--predict', action='store_true',
                        help="Generate prediction results with newly trained or input model.")
    parser.add_argument('--print_summary', action='store_true',
                        help="Display architecture of ML model.")
    parser.add_argument('--out_path', type=str, default=None,
                        help='Path to output reports with suggested format Cohort_Gene or Interval Name.')

    # Define variables from argument flags
    args = parser.parse_args()

    train_file = args.train_file
    test_file = args.test_file
    train_model = args.train
    feature_names = args.feature_names
    model_file = args.model_file
    predict = args.predict
    model_summary = args.print_summary
    out_path = args.out_path

    # Default feature names
    if len(feature_names) == 0:
        feature_names = ['dosage_interval', 'dosage_gene', 'del_dosage',
                         'std_baf', 'std_lrr', 'iqr_baf', 'iqr_lrr', 'avg_baf', 'avg_lrr']

    # ML data pre-processing for LSTM architecture
    X_train_reshaped, y_train, X_test_reshaped, train_samples, test_samples = prep_ml_datasets(
        feature_names, train_file, test_file)

    # Train model with reshaped training set
    if train_model:
        history = train_binary_lstm(X_train_reshaped, y_train, out_path)

    # Use newly trained or pre-trained model to make predictions on reshaped test sets
    if predict:
        if train_model:
            model_file = f'{out_path}_{X_train_reshaped.shape[1]}_windows.keras'
            print('Ready to predict using model:', model_file)
            model_predict(model_file, X_test_reshaped,
                          test_samples, out_path, model_summary)
        else:
            print('Ready to predict using model:', model_file)
            model_predict(model_file, X_test_reshaped,
                          test_samples, out_path, model_summary)


if __name__ == "__main__":
    main()
