import '../../styles/auth.css';

const AuthPage = () => {
  return (
    <div className="auth-page">
      <h1>Welcome to Omicsomics</h1>
      <p>Please sign in via your institutional identity provider.</p>
      <button type="button" className="auth-page__button">
        Sign in
      </button>
    </div>
  );
};

export default AuthPage;
