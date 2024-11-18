import pickle
import joblib
import pandas as pd
import numpy as np
import tensorflow as tf

# Supress Pandas copy warning
pd.options.mode.chained_assignment = None


def prep_ml_datasets(feature_names, train_path=None, test_path=None):
    """
    Prepares machine learning datasets for training and testing by reshaping them 
    into a 3D array format (samples, time steps per sample, number of features)
    suitable for an LSTM architecture.

    Arguments:
    feature_names (list of str): The names of the columns to be used as features.
    train_path (str, optional): Path to the CSV file containing the training data. Defaults to None.
    test_path (str, optional): Path to the CSV file containing the testing data. Defaults to None.

    Returns:
    tuple: A tuple containing:
        - X_train_reshaped (numpy.ndarray or None): The reshaped training feature array.
        - y_train (numpy.ndarray or None): The labels for the training data.
        - X_test_reshaped (numpy.ndarray or None): The reshaped testing feature array.
        - train_samples (list or None): List of sample IDs in the training set.
        - test_samples (list or None): List of sample IDs in the testing set.
    """
    if train_path:
        train_df = pd.read_csv(train_path)

        window_count = len(train_df.window.value_counts())
        print("Windows: ", window_count)

        X_train = train_df[feature_names]
        y_train = train_df[['IID', 'CNV_exists']][train_df['window'] == 0]
        y_train = y_train['CNV_exists'].values

        print("Training features:")
        print(X_train.shape)
        print("Training labels:")
        print(y_train.shape)

        # Extract unique sample IDs from the training set
        train_samples = list(train_df.IID[train_df['window'] == 0].values)

        # Reshape to 3D array (number of samples, time steps per sample, number of features)
        X_train_reshaped = X_train.to_numpy().reshape(
            (int(X_train.shape[0]/window_count), window_count, X_train.shape[1]))

        print("Reshaped training features:")
        print(X_train_reshaped.shape)
    else:
        y_train = None
        X_train_reshaped = None
        train_samples = None

    if test_path:
        test_df = pd.read_csv(test_path)

        window_count = len(test_df.window.value_counts())
        print("Windows: ", window_count)

        # Extract features for testing
        print("Testing features:")
        X_test = test_df[feature_names]
        print(X_test.shape)

        # Extract unique sample IDs from the testing set
        test_samples = list(test_df.IID[test_df['window'] == 0].values)
        X_test_reshaped = X_test.to_numpy().reshape(
            (int(X_test.shape[0]/window_count), window_count, X_test.shape[1]))

        print("Reshaped testing features:")
        print(X_test_reshaped.shape)
    else:
        X_test_reshaped = None
        test_samples = None

    return X_train_reshaped, y_train, X_test_reshaped, train_samples, test_samples


def train_binary_lstm(X_train_reshaped, y_train, out_path, verbosity=2, val_data=None):
    """
    Trains a binary LSTM (Long Short-Term Memory) model for classification tasks.
    The model outputs a single probability score per sample indicating the presence
    of specific CNV type (deletion or duplication). The trained model is saved in Keras format.

    Arguments:
    X_train_reshaped (numpy.ndarray): The reshaped 3D array for training features 
                                      with dimensions (samples, time steps, features).
    y_train (numpy.ndarray): The array of training labels.
    out_path (str): The path and filename (without extension) for saving the trained model.
    verbosity (int, optional): The verbosity level for training output (0, 1, or 2). Defaults to 2.
    val_data (tuple, optional): Validation data provided as (X_val, y_val) for evaluating 
                                the model during training. Defaults to None.

    Returns:
    keras.callbacks.History: The training history object containing metrics and loss over epochs.
    """
   # Build the binary sequence-to-vector LSTM model with 1 CNV type per model
    binary_lstm_model = tf.keras.models.Sequential([
        tf.keras.layers.LSTM(64, input_shape=(
            X_train_reshaped.shape[1], X_train_reshaped.shape[2]), return_sequences=True),
        tf.keras.layers.LSTM(64, return_sequences=True),
        tf.keras.layers.LSTM(64),
        tf.keras.layers.Dense(1, activation='hard_sigmoid')
    ])

    # Compile the model with appropriate loss function and metrics
    binary_lstm_model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy', 'AUC',
                                                                                     tf.keras.metrics.Precision(),
                                                                                     tf.keras.metrics.Recall()])

    # Fit model to reshaped training and validation data
    history = binary_lstm_model.fit(X_train_reshaped, y_train, batch_size=32, epochs=20,
                                    verbose=verbosity, validation_data=val_data)

    # Save the trained model in Keras format
    binary_lstm_model.save(f'{out_path}_windows.keras')

    # Save model - currently keras.src module issue in HPC
    # pickle.dump(binary_lstm_model, open(f'{out_path}_windows.sav', 'wb'))
    # joblib.dump(binary_lstm_model, f'{out_path}_windows.sav')

    return history


def model_predict(model_file, X_test_reshaped, test_samples, out_path, summary=True):
    """
    Loads a trained model, makes predictions on the test set, and saves the results to a CSV file.
    It also flags potential artifacts if more than 20% of the samples have a prediction probability
    over a certain threshold.

    Arguments:
    model_file (str): Path to the pre-trained model file (in Keras format).
    X_test_reshaped (numpy.ndarray): The reshaped 3D array for test features with 
                                     dimensions (samples, time steps, features).
    test_samples (list): List of sample IDs corresponding to the test data.
    out_path (str): Path and filename (without extension) for saving the prediction results.
    summary (bool, optional): Whether to print the model summary. Defaults to True.

    Returns:
    None: Outputs a CSV file containing the test results, including predicted 
          values and an artifact warning if applicable.
    """

    # Load in pre-trained model in Keras format
    loaded_model = tf.keras.models.load_model(model_file)

    # Load in model - currently keras.src module issue in HPC
    # loaded_model = pickle.load(file, 'rb')
    # loaded_model = pd.read_pickle(model_file)
    # loaded_model = joblib.load(open(model_file, 'rb'))

    # Print model summary if required
    if summary:
        print(loaded_model.summary())

    # Make predictions on the reshaped test set
    model_predictions = loaded_model.predict(X_test_reshaped)

    # Reshape predictions and create a DataFrame with results
    results_reshaped = model_predictions.reshape(-1)
    test_results = pd.DataFrame(
        {'IID': test_samples, 'Pred Values': results_reshaped})

    # Check for artifacts if over 20% of samples have a probability over 0.8
    test_results['AboveThreshold'] = test_results['Pred Values'] >= 0.8
    test_results['Artifact Warning'] = np.where(
        sum(test_results['AboveThreshold']) >= 0.2*len(test_results), 1, 0)
    test_results = test_results.drop(columns=['AboveThreshold'])

    # Print and save the results
    print(test_results)
    test_results.to_csv(f'{out_path}_windows_results.csv', index=False)
