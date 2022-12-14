"""DeepHF predictions"""
import numpy as np


def char2int(char):
    if char == 'A':
        return 1
    if char == 'T':
        return 2
    if char == 'C':
        return 3
    if char == 'G':
        return 4
    else:
        print('Received wrong char {} - exiting'.format(char))
        exit(1)


def get_predictions(model, config, sequences, decoded=False):
    if not decoded:
        sequences = prepare_arr(sequences)
    predictions = model.predict(sequences)
    # Scale outputs
    scaled_predictions = config.output_scale.inverse_transform(predictions)
    return scaled_predictions[:, 0]


def prepare_arr(sequences):
    decoded_seq = [decode_seq(seq) for seq in sequences]
    decoded_seq = np.concatenate(decoded_seq, axis=0)
    return decoded_seq


def decode_seq(seq):
    mer_array = np.array([0], dtype=np.uint8)
    for char in seq:
        num = char2int(char)
        mer_array = np.append(mer_array, np.array([num], dtype=np.uint8), axis=0)
    mer_array = np.expand_dims(mer_array, 0)
    return mer_array
