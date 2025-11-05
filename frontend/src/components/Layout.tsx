import { ReactNode } from 'react';
import Sidebar from './Sidebar';
import Topbar from './Topbar';
import '../styles/layout.css';

interface LayoutProps {
  children: ReactNode;
}

const Layout = ({ children }: LayoutProps) => {
  return (
    <div className="layout">
      <Sidebar />
      <div className="layout__main">
        <Topbar />
        <div className="layout__content">{children}</div>
      </div>
    </div>
  );
};

export default Layout;
