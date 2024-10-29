import pandas as pd
import numpy as np
import pickle
import joblib
import tensorflow as tf

from sklearn.preprocessing import StandardScaler, normalize

# Supress copy warning.
pd.options.mode.chained_assignment = None


# Eventually add CNN/other LSTM architectures
def prep_ml_datasets(feature_names, train_path=None, test_path=None):
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

        # Workaround for potentially multiple repeating samples
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

        # No testing labels - will check accuracy with manual observation for high probabiilties
        print("Testing features:")
        X_test = test_df[feature_names]
        print(X_test.shape)

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
    # Eventually add model that returns sequences & multiple probabilities for each CNV class
    binary_lstm_model = tf.keras.models.Sequential([
        tf.keras.layers.LSTM(64, input_shape=(
            X_train_reshaped.shape[1], X_train_reshaped.shape[2]), return_sequences=True),
        tf.keras.layers.LSTM(64, return_sequences=True),
        # not returning sequences b/c sequence-to-vector RNN model
        tf.keras.layers.LSTM(64),

        # dense final layer - right now: just outputting 1 class value per subject
        tf.keras.layers.Dense(1, activation='hard_sigmoid')
    ])

    # Will predict binary class value so need binary cross entropy for loss function
    binary_lstm_model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy', 'AUC',
                                                                                     tf.keras.metrics.Precision(),
                                                                                     tf.keras.metrics.Recall()])

    # Fit model to reshaped training vectors - make verbose option/more customizable architecture later
    history = binary_lstm_model.fit(X_train_reshaped, y_train, batch_size=32, epochs=20,
                                    verbose=verbosity, validation_data=val_data)  # default batch size = 32

    # Save model - keras.src module issue in swarm job
    # pickle.dump(binary_lstm_model, open(f'{out_path}_windows.sav', 'wb'))
    # joblib.dump(binary_lstm_model, f'{out_path}_windows.sav')

    binary_lstm_model.save(f'{out_path}_windows.keras')

    return history


def model_predict(model_file, X_test_reshaped, test_samples, out_path, summary=True):
    # Load in model - keras.src module issue in swarm job
    # loaded_model = pickle.load(file, 'rb')
    # loaded_model = pd.read_pickle(model_file)
    # loaded_model = joblib.load(open(model_file, 'rb'))

    # Load in model with keras file
    loaded_model = tf.keras.models.load_model(model_file)

    # If want summary
    if summary:
        print(loaded_model.summary())

    # Make predictions on reshaped test set
    model_predictions = loaded_model.predict(X_test_reshaped)

    # Display predicted classes and sample number to visually inspect
    results_reshaped = model_predictions.reshape(-1)
    test_results = pd.DataFrame(
        {'IID': test_samples, 'Pred Values': results_reshaped})

    # Binary value check for potential artifacts if over 20% of all samples have probability over 0.8
    test_results['AboveThreshold'] = test_results['Pred Values'] >= 0.8
    test_results['Artifact Warning'] = np.where(
        sum(test_results['AboveThreshold']) >= 0.2*len(test_results), 1, 0)
    test_results = test_results.drop(columns=['AboveThreshold'])

    print(test_results)
    test_results.to_csv(f'{out_path}_windows_results.csv', index=False)
