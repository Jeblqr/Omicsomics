import { useAuth } from '../contexts/AuthContext';

const Topbar = () => {
  const { user, logout } = useAuth();

  return (
    <header className="topbar">
      <div className="topbar__section">
        <span className="topbar__title">Unified Omics Platform</span>
      </div>
      <div className="topbar__section">
        <div className="topbar__user">
          {user?.full_name || user?.email || 'User'}
        </div>
        <button
          onClick={logout}
          className="logout-button"
          style={{
            marginLeft: '1rem',
            padding: '0.5rem 1rem',
            background: '#dc3545',
            color: 'white',
            border: 'none',
            borderRadius: '4px',
            cursor: 'pointer',
          }}
        >
          Logout
        </button>
      </div>
    </header>
  );
};

export default Topbar;
