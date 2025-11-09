from app.core import crypto


def main():
    user_id = 42
    msg = b"hello, this is a secret payload"
    blob = crypto.encrypt_for_user(user_id, msg)
    print("Encrypted blob len:", len(blob))
    pt = crypto.decrypt_for_user(user_id, blob)
    print("Decrypted matches:", pt == msg)


if __name__ == "__main__":
    main()
