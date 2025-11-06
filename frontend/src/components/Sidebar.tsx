import { NavLink } from 'react-router-dom';

const Sidebar = () => {
  const links = [
    { to: '/', label: 'Dashboard' },
    { to: '/projects', label: 'Projects' },
    { to: '/runs', label: 'Runs' },
    { to: '/data', label: 'Data' },
    { to: '/settings', label: 'Settings' },
  ];

  return (
    <aside className="sidebar">
      <h1 className="sidebar__logo">Omicsomics</h1>
      <nav className="sidebar__nav">
        {links.map((link) => (
          <NavLink
            key={link.to}
            to={link.to}
            className={({ isActive }) =>
              `sidebar__link${isActive ? ' sidebar__link--active' : ''}`
            }
            end={link.to === '/'}
          >
            {link.label}
          </NavLink>
        ))}
      </nav>
    </aside>
  );
};

export default Sidebar;
